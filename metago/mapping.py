import sys
import os
import random
import string

import click
from pyflux import FluxWorkflowRunner
from subprocess import call


class SampleReadPileup(FluxWorkflowRunner):
    """ Coordinates sample read mapping and pileup visualization creation """
    def __init__(self, fastq, reference, visualize, output, max_ppn, max_mem, paired):
        self.fastq = fastq
        self.reference = reference
        self.visualize = visualize
        self.output = output
        self.max_ppn = max_ppn
        self.max_mem = max_mem
        self.paired = paired

    def workflow(self):
        fp = os.path.dirname(os.path.abspath(__file__))
        env = os.path.join(fp, '../dependencies/miniconda/bin/activate')
        output_dp = os.path.join(self.output, 'read_pileup')
        if not os.path.exists(output_dp): os.makedirs(output_dp)

        submitted_cmds = []

        sam_fp = os.path.join(output_dp, '{}.sam'.format(self.reference.split('/')[-1].split('.')[0]))        
        if not os.path.exists(sam_fp):
            cmd = 'source {} && '.format(env) 
            cmd += 'bbmap.sh -Xmx{}m in={} ref={} t={} outm={}'.format(self.fastq, self.max_mem, self.reference, self.max_ppn, sam_fp)
            if not self.paired:
                cmd += ' interleaved=f'
            self.addTask('bbmap', nCores=self.max_ppn, memMb=self.max_mem, command=cmd)
            submitted_cmds.append('bbmap')

        bam_fp = sam_fp.replace('.sam', '.bam')
        if not os.path.exists(bam_fp):
            cmd = 'source {} && '.format(env)
            cmd += 'samtools sort {} > {}'.format(sam_fp, bam_fp)
            print '\n\n{}\n\n'.format(cmd)
            # based on samtools sort default memory per thread
            self.addTask('samsort', nCores=1, memMb=768, command=cmd, dependencies=submitted_cmds)
            submitted_cmds.append('samsort')

        index_fp = bam_fp + '.bai'
        if not os.path.exists(index_fp):
            cmd = 'source {} && '.format(env)
            cmd += 'samtools index {}'.format(bam_fp)
            self.addTask('samindex', nCores=1, memMb=768, command=cmd, dependencies=submitted_cmds)
            submitted_cmds.append('samindex')

        if self.visualize:
            cmd = 'source {} && '.format(env)
            cmd += 'pyleup visualize {} -o {}'.format(bam_fp, output_dp) 
            self.addTask('pileup', nCores=1, memMb=768, command=cmd, dependencies=submitted_cmds)
        


class RunReadPileup(FluxWorkflowRunner):
    """ Coordinate multiple sample pileups from a run directory
    
    Assumes samples have already been processed by quality control step.
    """
    def __init__(self, run_dp, reference, visualize, output, max_ppn, max_mem, paired):
        self.run_dp = run_dp
        self.reference = reference
        self.visualize = visualize
        self.output = output
        self.max_ppn = max_ppn
        self.max_mem = max_mem
        self.paired = paired

    def workflow(self):
        for sample in os.listdir(self.run_dp):
            sample_dp = os.path.join(self.run_dp, sample)
            if not os.path.isdir(sample_dp) or not 'Sample_' in sample: continue
            sample_fastq = None
            for fastq in os.listdir(sample_dp):
                if '.fastq' not in fastq and not '.fa' in fastq: continue
                sample_fastq = os.path.join(sample_dp, fastq)
                break
            sample_output = os.path.join(self.output, sample)
            sample_runner = SampleReadPileup(fastq=sample_fastq,
                                             reference=self.reference,
                                             visualize=self.visualize,
                                             output=sample_output,
                                             max_ppn=self.max_ppn,
                                             max_mem=self.max_mem,
                                             paired=self.paired)
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
@click.option('--run', help='full path to ruin directory with quality controlled samples')
@click.option('--sample', help='full path to fastq')
@click.option('--reference', '-r', required=True)
@click.option('--visualize/--not-visualize', default=False)
@click.option('--paired/--not-paired', default=False)
@click.pass_context
def mapping(ctx, run, sample, reference, visualize, paired):
    """ Coordinate mapping of reads to reference """
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output = os.path.join(ctx.obj['OUTPUT'], 'logs/pileup_{}'.format(r))

    if run:
        runner = RunReadPileup(run_dp=run,
                           reference=reference,
                           visualize=visualize,
                           output=ctx.obj['OUTPUT'],
                           max_ppn=ctx.obj['PPN'],
                           max_mem=ctx.obj['MEM'],
                           paired=paired)
    elif sample:
        runner = SampleReadPileup(fastq=sample,
                           reference=reference,
                           visualize=visualize,
                           output=ctx.obj['OUTPUT'],
                           max_ppn=ctx.obj['PPN'],
                           max_mem=ctx.obj['MEM'],
                           paired=paired)
    return_code = runner.run(mode='local', dataDirRoot=log_output, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])
    if return_code != 0: sys.exit('Non-Zero Exit Code in run_qc')
    return return_code


if __name__ == "__main__":
    cli(obj={})
        
