import sys
import os
import random
import string

import click
from pyflux import FluxWorkflowRunner
from subprocess import call


class SampleBinning(FluxWorkflowRunner):
    """ Follows Metahit Assembly Tips

    https://github.com/voutcn/megahit/wiki/Assembly-Tips
    """
    def __init__(self, contig_fp, fastq_fp, output_dp, max_ppn, max_mem):
        self.contig_fp = contig_fp
        self.fastq_fp = fastq_fp
        self.output_dp = os.path.join(output_dp, 'binning')

        self.max_ppn = int(max_ppn)
        self.max_mem = int(max_mem)


    def workflow(self):
        """ Sample Assembly Workflow
        
        To consider: estimating diversity (genome complexity) and depth to inform assembly parameter settings
        """
        fp = os.path.dirname(os.path.abspath(__file__))
        conda = os.path.join(fp, '../dependencies/miniconda/bin/activate')

        sample_output_dp = os.path.join(self.output_dp, self.sid)
        if not os.path.exists(sample_output_dp): os.makedirs(sample_output_dp)

        if not os.path.exists(assembly_dp):
            os.makedirs(assembly_dp)
        if len(os.listdir(assembly_dp)) == 0:
            cmd = 'source {} && run_MaxBin.pl -contig {} -reads {} -thread {} -out {}'.format(
                             conda, self.contig_fp, self.fastq_fp, self.max_ppn, sample_output_dp)
            print 'cmd: {}'.format(cmd)
            self.addTask('bin', nCores=self.max_ppn, memMb=self.max_mem, command=cmd)


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
@click.option('--contig_fp', '-c', required=True, help='file path to contig fasta file')
@click.option('--read_fp', '-r', required=True, help='file path to read fasta/fastq file')
@click.pass_context
def sample_binning(ctx, sample_dp):
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(ctx.obj['OUTPUT'], 'logs/sample_assembly_{}'.format(r))

    sample_fastq = ''
    for fastq in os.listdir(sample_dp):
        if 'fastq' not in fastq: continue
        fastq_fp = os.path.join(sample_dp, fastq)
        sample_fastq = fastq_fp
        break # this file should be a single interleaved and quality controlled fastq

    runner = SampleBinning(contig_fp=contig_fp, fastq_fp=read_fp, output_dp=ctx.obj['OUTPUT'],
                                                max_ppn=ctx.obj['PPN'], max_mem=ctx.obj['MEM'])
    runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM']) 


if __name__ == "__main__":
    cli(obj={})
