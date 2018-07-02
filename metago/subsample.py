import sys
import os
import random
import string

import click
from pyflux import FluxWorkflowRunner
from subprocess import call


class SampleSubSampler(FluxWorkflowRunner):
    def __init__(self, fasta, min_ratio, max_ratio, ratio_interval, output_dp, max_ppn, max_mem):
        self.fasta = fasta
        self.min_ratio = min_ratio
        self.max_ratio = max_ratio
        self.ratio_interval = ratio_interval
        self.output_dp = os.path.join(output_dp, 'subsampling_output')

        self.max_ppn = int(max_ppn)
        self.max_mem = int(max_mem)


    def workflow(self):
        fp = os.path.dirname(os.path.abspath(__file__))
        conda = os.path.join(fp, '../dependencies/miniconda/bin/activate')

        for x in np.arange(self.min_ratio, self.max_ratio+self.ratio_interval, self.ratio_interval):
            cmd = 'source {} && reformat.sh in={} out={}/{}p_{} samplerate={}'.format(
                          conda, self.fasta, self.output_dp, int(x*10), os.path.basename(self.fasta))
            print 'cmd: {}'.format(cmd)
            self.addTask('subsample_{}'.format(x), nCores=self.max_ppn, memMb=self.max_mem, command=cmd)


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


@click.pass_context
@click.argument('fasta_fp')
@click.option('--min_ratio', default=0.1)
@click.option('--max_ratio', default=0.9)
@click.option('--ratio_interval', default=0.1)
def subsample(ctx, fasta_fp):
    """ VirSorter subworkflow manager

    Arguments:
    fasta_fp -- file path to contig file path to be used for viral sorting
    """
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(ctx.obj['OUTPUT'], 'logs/minhash_{}'.format(r))

    runner = SampleSubSampler(fasta=fasta_fp, min_ratio=min_ratio, max_ratio=max_ratio, ratio_interval=ratio_interval,
                                        output_dp=ctx.obj['OUTPUT'], max_ppn=ctx.obj['PPN'], max_mem=ctx.obj['MEM'])
    runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])


if __name__ == "__main__":
    cli(obj={})
