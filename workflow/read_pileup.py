import sys
import os
import random
import string

import click
from pyflux import FluxWorkflowRunner
from subprocess import call


class SampleReadPileup(FluxWorkflowRunner):
    """ Coordinates sample read mapping and pileup visualization creation """
    def __init__(self, fastq, referece, output, max_ppn, max_mem):
        self.fastq = fastq
        self.reference = reference
        self.output = output
        self.max_ppn = max_ppn
        self.max_mem = max_mem

    def workflow(self):
        fp = os.path.dirname(os.path.abspath(__file__))
        env = os.path.join(fp, '../dependencies/miniconda/bin/activate')
        output_dp - os.path.join(self.output, 'read_pileup')
        if not os.path.exists(output_dp): os.makedirs(output_dp)
        
        cmd = 'source {} && '.format(env) 
        cmd += 'bbmap in={} ref={} t={}'.format(self.fastq, self.reference, self.max_ppn)
        self.addTask('bbmap', nCores=self.max_ppn, memMb=self.max_mem)


class RunReadPileup(FluxWorkflowRunner):
    """ Coordinate multiple sample pileups from a run directory
    
    Assumes samples have already been processed by quality control step.
    """
    def __init__(self, run_dp, reference, output, max_ppn, max_mem):
        self.run_dp = run_dp
        self.reference = reference
        self.output = output
        self.max_ppn = max_ppn
        self.max_mem = max_mem

    def workflow(self):
        for sample in os.listdir(self.run_dp):
            sample_dp = os.path.join(self.run_dp, sample)
            if not os.path.isdir(sample_dp) or not 'Sample_' in sample: continue
            sample_fastq = None
            for fastq in os.listdir(sample_dp):
                if '.fastq' not in fastq and not '.fa' in fastq: continue
                sample_fastq = os.path.join(sample_dp, fastq)
                break
            sample_runner = SampleReadPileup(fastq=sample_fastq, output=self.output, max_ppn=self.max_ppn, max_mem=self.max_mem)
            self.addWorkflowTask(label=sample, workflowRunnerInstance=sample_runner)

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
@click.option('--reference', '-r', required=True)
@click.pass_context
def run_pileup(ctx, run_dp, reference):
    """ Coordinate mapping of reads from run samples to reference and visualize """
    pass


@cli.command()
@click.argument('sample_dp')
@click.option('--reference', '-r', required=True)
@click.pass_context
def sample_pileup(ctx, sample_dp):
    """ Map reads from sample fastqs/fastas to reference and visualize """
    pass


@cli.command()
@click.argument('read_fp')
@click.option('--reference', '-r', required=True)
@click.pass_context
def read_pileup(ctx, read_fp, reference):
    """ Map reads from fastq/fasta to reference and visualize """
    pass

if __name__ == "__main__":
    cli(obj={})
        
