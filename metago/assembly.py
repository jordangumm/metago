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
    def __init__(self, sid, fastq, output_dp, max_ppn, max_mem):
        self.sid = sid
        self.fastq = fastq
        self.output_dp = output_dp

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

        scheduled_tasks = []
        normalized_fp = os.path.join(sample_output_dp, 'normalized', '{}.fastq'.format(self.sid))
        if not os.path.exists(normalized_fp):
            cmd = 'source {} && bbnorm.sh -Xmx{}m t={}'.format(conda, self.max_mem-2000, self.max_ppn-1)
            cmd += ' in={} out={} target=100 min=5'.format(self.fastq, normalized_fp)
            print 'cmd: {}'.format(cmd)
            self.addTask("normalize_{}".format(self.sid), nCores=self.max_ppn, memMb=self.max_mem, command=cmd)
            scheduled_tasks.append("normalize_{}".format(self.sid))

        assembly_dp = os.path.join(sample_output_dp, 'assembly')
        if not os.path.exists(assembly_dp):
            os.makedirs(assembly_dp)
        if len(os.listdir(assembly_dp)) == 0:
            cmd = 'source {} && megahit --12 {} -t {} -o {}'.format(conda, normalized_fp, self.max_cpu, assembly_dp)
            print 'cmd: {}'.format(cmd)
            self.addTask('assemble_{}'.format(self.sid), nCores=self.max_ppn, memMb=self.max_mem, command=cmd)
            scheduled_tasks.append('assemble_{}'.format(self.sid))


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
@click.argument('sample_dp')
@click.pass_context
def sample_assembly(ctx, sample_dp):
    """ Sample assembly subworkflow manager """
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(ctx.obj['OUTPUT'], 'logs/sample_assembly_{}'.format(r))

    sid = os.path.basename(sample_dp).replace('_Sample','')

    sample_fastq = ''
    for fastq in os.listdir(sample_dp):
        if 'fastq' not in fastq: continue
        fastq_fp = os.path.join(sample_dp, fastq)
        sample_fastq = fastq_fp
        break # this file should be a single interleaved and quality controlled fastq

    runner = SampleAssembly(sid=sid, fastq=sample_fastq, output_dp=ctx.obj['OUTPUT'],
                            max_ppn=ctx.obj['PPN'], max_mem=ctx.obj['MEM'])
    runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM']) 


@cli.command()
@click.argument('run_dp')
@click.pass_context
def run_assembly(ctx, run_dp):
    """ Run assembly subworkflow manager

    Arguments:
    run_dp -- String path to run directory to use for analysis
    """
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(ctx.obj['OUTPUT'], 'logs/run_assembly_{}'.format(r))

    runner = RunSampleAssembly(run_dp=run_dp, output_dp=ctx.obj['OUTPUT'],
                               max_ppn=ctx.obj['PPN'], max_mem=ctx.obj['MEM'])
    runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])


if __name__ == "__main__":
    cli(obj={})
