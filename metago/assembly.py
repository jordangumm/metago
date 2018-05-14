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
        """ Quality Control Workflow
        
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

            sample_assembly_runner = SampleAssembly(sid=sample.replace('Sample_',''), fastq=sample_fastq, output_dp=self.output_dp,
                                                                     max_ppn=self.max_ppn, max_mem=self.max_mem)
            self.addWorkflowTask(label=sample, workflowRunnerInstance=sample_assembly_runner)


@click.command()
@click.option('--run_dp', '-r', required=True)
@click.option('--output_dp', '-o', required=True)
@click.option('--ppn', '-p', required=True)
@click.option('--mem', '-m', required=True)
def runner(run_dp, output_dp, ppn, mem):
    """ Analysis Workflow Management

    Sets up Pyflow WorkflowRunner and launches locally by default or via flux

    Arguments:
    run_dp -- String path to run directory to use for analysis
    """
    log_output_dp = os.path.join(output_dp, 'logs')

    workflow_runner = RunSampleAssembly(run_dp=run_dp, output_dp=output_dp, max_ppn=ppn, max_mem=mem)
    workflow_runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ppn, memMb=mem)


if __name__ == "__main__":
    runner()
