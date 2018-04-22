#!/usr/bin/env python

import os, sys
import click

from workflow.pyflux import FluxWorkflowRunner
from workflow.quality_control import RunQualityControl
from subprocess import call


class Runner(FluxWorkflowRunner):
    def __init__(self, run_dp, output_dp, max_ppn, max_mem):
        self.run_dp = run_dp
        self.output_dp = output_dp
        self.max_ppn = max_ppn
        self.max_mem = max_mem

    def workflow(self):
        """ method invoked on class instance run call """
        #qc_runner = RunQualityControl(run_dp=self.run_dp, output_dp=self.output_dp, max_ppn=self.max_ppn, max_mem=self.max_mem)
        fp = os.path.dirname(os.path.abspath(__file__))
        conda = os.path.join(fp, 'dependencies/miniconda/bin/activate')
        qc_runner = os.path.join(fp, 'workflow/quality_control.py')
        self.addTask("qc", nCores=self.max_ppn, memMb=self.max_mem,
                     command='source {} && python {} -r {} -o {} -p {} -m {}'.format(conda,
                                qc_runner, self.run_dp, self.output_dp, self.max_ppn, self.max_mem))
        
        #self.addWorkflowTask(label="qc", workflowRunnerInstance=qc_runner)


@click.command()
@click.argument('run_dp')
@click.option('--output', '-o', default=None)
@click.option('--flux/--no-flux', default=False)
@click.option('--dispatch/--no-dispatch', default=True)
@click.option('--account', '-a')
@click.option('--ppn', '-p', default=4)
@click.option('--mem', '-m', default='20000') # current limitation, only handles mb
@click.option('--walltime', '-w', default='2:00:00')
def runner(run_dp, output, flux, dispatch, account, ppn, mem, walltime):
    """ Analysis Workflow Management

    Sets up Pyflow WorkflowRunner and launches locally by default or via flux

    Arguments:
    run_dp -- String path to run directory of Illumina Sample directories to use in analysis
    """
    if not output: output = os.path.join(run_dp, 'bioinfo')
    log_output_dp = os.path.join(output, 'logs', 'runner')
    workflow_runner = Runner(run_dp=run_dp, output_dp=output, max_ppn=ppn, max_mem=mem)

    if flux:
        if not account: sys.exit('To attempt a submission to the flux cluster you need to supply an --account/-a')
        #if dispatch:
        full_dp = os.path.dirname(os.path.abspath(__file__))
        activate = 'source {}'.format(os.path.join(full_dp, 'dependencies', 'miniconda', 'bin', 'activate'))
        runner_fp = os.path.join(full_dp, 'runner.py')
        qsub = 'qsub -N pyflux_handler -A {} -q fluxm -l nodes=1:ppn={},mem={}mb,walltime={}'.format(
                                                                   account, ppn, mem, walltime)
        call('echo "{} && python {} {} -o {} -a {} -p {} -m {} -w {}" | {}'.format(
                                                                   activate, runner_fp,
                                                                   run_dp, output, account,
                                                                   ppn, mem, walltime, qsub), shell=True)
       # else:
       #     workflow_runner.run(mode='flux', dataDirRoot=log_output_dp, nCores=ppn, memMb=mem,
       #                         schedulerArgList=['-N', 'viral_runner',
       #                                           '-A', account,
       #                                           '-l', 'nodes=1:ppn={},mem={}mb,walltime={}'.format(ppn, mem, walltime)])
    else:
        workflow_runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ppn, memMb=mem)


if __name__ == "__main__":
    runner()
