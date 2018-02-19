import sys
import os
import random
import string

import click
from pyflow import WorkflowRunner
from subprocess import call

class Runner(WorkflowRunner):
    def __init__(self, run_dp, num_cpu):
        self.run_dp = run_dp

    def workflow(self):
        """ method invoked on class instance run call """
        self.addTask("tmp", command=['echo', '"Viral Omics Analysis Workflow: Running Under {}"'.format(self.run_dp)])


@click.command()
@click.argument('run_dp')
@click.option('--flux/--no-flux', default=False)
@click.option('--account', '-a')
@click.option('--ppn', '-p', default=8)
@click.option('--mem', '-m', default='20gb')
@click.option('--walltime', '-w', default='2:00:00')
def runner(run_dp, flux, account, ppn, mem, walltime):
    """ Analysis Workflow Management

    Sets up Pyflow WorkflowRunner and launches locally by default or via flux

    Arguments:
    run_dp -- String path to run directory to use for analysis
    """
    log_output_dp = os.path.join(run_dp, 'bioinfo', 'logs', 'runner')

    if flux:
        full_dp = os.path.dirname(os.path.abspath(__file__))
        activate = 'source {}'.format(os.path.join(full_dp, 'dependencies', 'miniconda', 'bin', 'activate'))
        runner_fp = os.path.join(full_dp, 'runner.py')
        qsub = 'qsub -N omics_16s -A {} -q fluxm -l nodes=1:ppn={},mem={},walltime={}'.format(account, ppn, mem, walltime)
        call('echo "{} && python {} {}" | {}'.format(activate, runner_fp, run_dp, qsub), shell=True)
    else:
        workflow_runner = Runner(run_dp=run_dp, num_cpu=ppn)
        workflow_runner.run(mode='local', dataDirRoot=log_output_dp)


if __name__ == "__main__":
    runner()
