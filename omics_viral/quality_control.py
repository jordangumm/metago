import sys
import os
import random
import string

import click
from pyflow import WorkflowRunner
from subprocess import call

class QualityControl(WorkflowRunner):
    def __init__(self, run_dp, num_cpu):
        self.run_dp = run_dp

    def workflow(self):
        """ method invoked on class instance run call """
        self.addTask("tmp", command=['echo', '"Viral Omics Analysis Workflow: QC Under {}"'.format(self.run_dp)])


def runner(run_dp, flux, account, ppn, mem, walltime):
    """ Analysis Workflow Management

    Sets up Pyflow WorkflowRunner and launches locally by default or via flux

    Arguments:
    run_dp -- String path to run directory to use for analysis
    """
    log_output_dp = os.path.join(run_dp, 'bioinfo', 'logs', 'quality_control')

    workflow_runner = Runner(run_dp=run_dp, num_cpu=ppn)
    workflow_runner.run(mode='local', dataDirRoot=log_output_dp)


if __name__ == "__main__":
    runner()
