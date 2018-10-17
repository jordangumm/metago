import sys
import os
import random
import string

import click
from pyflux import FluxWorkflowRunner
from subprocess import call



class CoAssembly(FluxWorkflowRunner):
    pass


class SampleAssembly(FluxWorkflowRunner):
    """ Follows Metahit Assembly Tips

    https://github.com/voutcn/megahit/wiki/Assembly-Tips
    """
    def __init__(self, sid, fastq, output_dp, max_ppn, max_mem, continue_assembly=False):
        self.sid = sid
        self.fastq = fastq
        self.output_dp = output_dp

        self.max_ppn = int(max_ppn)
        self.max_mem = int(max_mem)

        self.continue_assembly = continue_assembly


    def workflow(self):
        """ Sample Assembly Workflow
        
        To consider: estimating diversity (genome complexity) and depth to inform assembly parameter settings
        """
        fp = os.path.dirname(os.path.abspath(__file__))
        conda = os.path.join(fp, '../dependencies/miniconda/bin/activate')

        sample_output_dp = os.path.join(self.output_dp, self.sid)
        if not os.path.exists(sample_output_dp): os.makedirs(sample_output_dp)

        assembly_dp = os.path.join(sample_output_dp, 'assembly')
        if not os.path.exists(assembly_dp):
            cmd = 'source {} && megahit --preset meta-sensitive --12 {} -t {} -o {}'.format(conda, self.fastq, self.max_ppn, assembly_dp)
            print 'cmd: {}'.format(cmd)
            self.addTask('assemble_{}'.format(self.sid), nCores=self.max_ppn, memMb=self.max_mem, command=cmd)
        elif self.continue_assembly:
            cmd = 'source {} && megahit --preset meta-sensitive --12 {} -t {} -o {} --continue'.format(
                                                                     conda, self.fastq, self.max_ppn, assembly_dp)
            print 'cmd: {}'.format(cmd)
            self.addTask('assemble_{}'.format(self.sid), nCores=self.max_ppn, memMb=self.max_mem, command=cmd)


class RunSampleAssembly(FluxWorkflowRunner):
    def __init__(self, run_dp, output_dp, max_ppn, max_mem):
        self.run_dp = run_dp
        self.output_dp = output_dp

        self.max_ppn = max_ppn
        self.max_mem = max_mem

    def workflow(self):
        for sample in os.listdir(self.run_dp):
            sample_dp = os.path.join(self.run_dp, sample)
            if not os.path.isdir(sample_dp) or not 'Sample_' in sample: continue
            sample_fastq = ''
            for fastq in os.listdir(sample_dp):
                if 'fastq' not in fastq: continue
                fastq_fp = os.path.join(sample_dp, fastq)
                sample_fastq = fastq_fp
                break # this file should be a single interleaved and quality controlled fastq

            sample_assembly_runner = SampleAssembly(sid=sample.replace('Sample_',''),
                                                    fastq=sample_fastq,
                                                    output_dp=self.output_dp,
                                                    max_ppn=self.max_ppn,
                                                    max_mem=self.max_mem)
            self.addWorkflowTask(label=sample, workflowRunnerInstance=sample_assembly_runner)


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
@click.option('--fastqs', '-f', multiple=True)
@click.option('--continue_assembly/--no_continue_assembly', default=False)
@click.pass_context
def co_assembly(ctx, fastqs, continue_assembly):
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(ctx.obj['OUTPUT'], 'logs/sample_assembly_{}'.format(r))

    if len(fastqs) < 2:
        sys.exit('[ERROR]: Not enough fastq files entered')
    fastqs = ','.join(fastqs)

    runner = SampleAssembly(sid='co_assembly', fastq=fastqs, output_dp=ctx.obj['OUTPUT'],
                        max_ppn=ctx.obj['PPN'], max_mem=ctx.obj['MEM'], continue_assembly=continue_assembly)
    runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])
    
    


@cli.command()
@click.option('--sample', help='directory path to sample that has been quality controlled')
@click.option('--run', help='directory path to run with samples that have been quality controlled')
@click.option('--continue_assembly/--no_continue_assembly', default=False)
@click.pass_context
def assembly(ctx, sample, run, continue_assembly):
    """ Assembly subworkflow manager """
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(ctx.obj['OUTPUT'], 'logs/sample_assembly_{}'.format(r))

    if sample:
        sid = os.path.basename(sample).replace('_Sample','')

        sample_fastq = ''
        for fastq in os.listdir(sample):
            if 'fastq' not in fastq: continue
            fastq_fp = os.path.join(sample, fastq)
            sample_fastq = fastq_fp
            break # this file should be a single interleaved and quality controlled fastq

        runner = SampleAssembly(sid=sid, fastq=sample_fastq, output_dp=ctx.obj['OUTPUT'],
                            max_ppn=ctx.obj['PPN'], max_mem=ctx.obj['MEM'], continue_assembly=continue_assembly)
        runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM']) 
    elif run:
        runner = RunSampleAssembly(run_dp=run, output_dp=ctx.obj['OUTPUT'],
                               max_ppn=ctx.obj['PPN'], max_mem=ctx.obj['MEM'])
        runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])


if __name__ == "__main__":
    cli(obj={})
