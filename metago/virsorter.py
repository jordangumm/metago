import sys
import os
import random
import string

import click
from pyflux import FluxWorkflowRunner
from subprocess import call


class VirSorter(FluxWorkflowRunner):
    def __init__(self, fasta, output_dp, max_ppn, max_mem):
        self.fasta = fasta
        self.output_dp = os.path.join(output_dp, 'virsorter_output')

        self.max_ppn = int(max_ppn)
        self.max_mem = int(max_mem)


    def workflow(self):
        fp = os.path.dirname(os.path.abspath(__file__))
        conda = os.path.join(fp, '../dependencies/miniconda/bin/activate')
        data_dir = os.path.join(fp, '../dependencies/virsorter-data')

        cmd = 'source {} && wrapper_phage_contigs_sorter_iPlant.pl --no_c --virome -f {} --db 1 --wdir {} --ncpu {} --data-dir {}'.format(
                                                                 conda, self.fasta, self.output_dp, self.max_ppn, data_dir)
        print 'cmd: {}'.format(cmd)
        self.addTask('virsorter', nCores=self.max_ppn, memMb=self.max_mem, command=cmd)


@click.group()
@click.option('--output', '-o', required=True)
@click.option('--ppn', '-p', required=True)
@click.option('--mem', '-m', required=True)
@click.pass_context
def cli(ctx, output, ppn, mem):
    if not os.path.exists(output): os.makedirs(output)
    ctx.obj['OUTPUT'] = output
    ctx.obj['PPN'] = ppn
    ctx.obj['MEM'] = mem


@cli.command()
@click.pass_context
@click.argument('fasta_fp')
def virsorter(ctx, fasta_fp):
    """ VirSorter subworkflow manager

    Arguments:
    fasta_fp -- file path to contig file path to be used for viral sorting
    """
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(ctx.obj['OUTPUT'], 'logs/minhash_{}'.format(r))

    runner = VirSorter(fasta=fasta_fp, output_dp=ctx.obj['OUTPUT'], max_ppn=ctx.obj['PPN'], max_mem=ctx.obj['MEM'])
    runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])


if __name__ == "__main__":
    cli(obj={})
