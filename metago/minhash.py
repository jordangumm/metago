import sys
import os
import random
import string

import click
from pyflux import FluxWorkflowRunner
from subprocess import call




class SampleCompute(FluxWorkflowRunner):
    def __init__(self, fastq, output_dp, num_hashes, kmer_length, max_ppn, max_mem):
        self.fastq = fastq
        self.output_dp = output_dp
        self.num_hashes = num_hashes
        self.kmer_length = kmer_length

        self.max_ppn = int(max_ppn)
        self.max_mem = int(max_mem)


    def workflow(self):
        """ Sample Assembly Workflow
        
        To consider: estimating diversity (genome complexity) and depth to inform assembly parameter settings
        """
        fp = os.path.dirname(os.path.abspath(__file__))
        conda = os.path.join(fp, '../dependencies/miniconda/bin/activate')

        output_fp = '{}.sig'.format(os.path.join(self.output_dp, self.fastq.split('/')[-1].split('.')[0]))
        if not os.path.exists(output_fp):
            cmd = 'source {} && source activate py3 && sourmash compute {} --scaled {} -k {} -o {}'.format(
                                            conda, self.fastq, self.num_hashes, self.kmer_length, output_fp)
            print 'cmd: {}'.format(cmd)
            self.addTask("compute", nCores=1, memMb=1000, command=cmd)


class RunMinHash(FluxWorkflowRunner):
    def __init__(self, run_dp, output_dp, num_hashes, kmer_length, max_ppn, max_mem):
        self.run_dp = run_dp
        self.output_dp = output_dp
        self.num_hashes = num_hashes
        self.kmer_length = kmer_length

        self.max_ppn = max_ppn
        self.max_mem = max_mem

    def workflow(self):
        fp = os.path.dirname(os.path.abspath(__file__))
        conda = os.path.join(fp, '../dependencies/miniconda/bin/activate')

        scheduled_jobs = []

        signatures_dp = os.path.join(self.output_dp, 'sourmash_signatures')
        if not os.path.exists(signatures_dp): os.makedirs(signatures_dp)
        for sample in os.listdir(self.run_dp):
            sample_dp = os.path.join(self.run_dp, sample)
            if not os.path.isdir(sample_dp) or not 'Sample_' in sample: continue
            sample_fastq = ''
            for fastq in os.listdir(sample_dp):
                if 'fastq' not in fastq: continue
                fastq_fp = os.path.join(sample_dp, fastq)
                sample_fastq = fastq_fp
                break # this file should be a single interleaved and quality controlled fastq

            sample_assembly_runner = SampleCompute(fastq=sample_fastq,
                                                    output_dp=signatures_dp,
                                                    num_hashes=self.num_hashes,
                                                    kmer_length=self.kmer_length,
                                                    max_ppn=self.max_ppn,
                                                    max_mem=self.max_mem)
            self.addWorkflowTask(label=sample, workflowRunnerInstance=sample_assembly_runner)
            scheduled_jobs.append(sample)

        compare_fp = os.path.join(self.output_dp, 'sourmash_compare')
        if not os.path.exists(compare_fp):
            cmd = 'source {} && source activate py3 && cd {} && sourmash compare *.sig -o {}'.format(
                                                                          conda, signatures_dp, compare_fp)
            self.addTask('compare', nCores=self.max_ppn, memMb=self.max_mem, command=cmd, dependencies=scheduled_jobs)
            scheduled_jobs.append('compare')


        compare_csv_fp = os.path.join(self.output_dp, 'sourmash_compare.csv')
        if not os.path.exists(compare_fp):
            cmd = 'source {} && source activate py3 && cd {} && sourmash compare *.sig --csv {}'.format(
                                                                          conda, signatures_dp, compare_csv_fp)
            self.addTask('compare_csv', nCores=self.max_ppn, memMb=self.max_mem, command=cmd, dependencies=scheduled_jobs)
            scheduled_jobs.append('compare_csv')
        

        plot_dp = os.path.join(self.output_dp, 'sourmash_plot')
        cmd = 'source {} && source activate py3 && cd {} && sourmash plot {}'.format(conda, plot_dp, compare_fp)
        self.addTask('plot', nCores=self.max_ppn, memMb=self.max_mem, command=cmd, dependencies=scheduled_jobs)


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
@click.argument('run_dp')
@click.option('--num_hashes', '-s', default=1000)
@click.option('--kmer_length', '-k', default=31, help='minhash kmer length to use (default=31, recommended for genus level)')
@click.pass_context
def run_minhash(ctx, run_dp, num_hashes, kmer_length):
    """ Run assembly subworkflow manager

    Recommended number of hashes is --scaled 1000, which is the default here.
    http://sourmash.readthedocs.io/en/latest/using-sourmash-a-guide.html#what-resolution-should-my-signatures-be-how-should-i-compute-them

    Arguments:
    run_dp -- String path to run directory to use for analysis
    """
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(ctx.obj['OUTPUT'], 'logs/minhash_{}'.format(r))

    runner = RunMinHash(run_dp=run_dp, output_dp=ctx.obj['OUTPUT'], num_hashes=num_hashes, kmer_length=kmer_length,
                                                           max_ppn=ctx.obj['PPN'], max_mem=ctx.obj['MEM'])
    runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])


if __name__ == "__main__":
    cli(obj={})
