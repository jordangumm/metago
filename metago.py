#!/usr/bin/env python

import os, sys
import click
import random
import string

from workflow.pyflux import FluxWorkflowRunner
from workflow.quality_control import RunQualityControl
from workflow.quality_control import SampleQualityControl
from subprocess import call


class Runner(FluxWorkflowRunner):
    def __init__(self, cmd, output_dp, max_ppn, max_mem):
        self.cmd = cmd
        self.output_dp = output_dp
        self.max_ppn = max_ppn
        self.max_mem = max_mem

    def workflow(self):
        """ method invoked on class instance run call """
        self.addTask("runner", nCores=self.max_ppn, memMb=self.max_mem, command=self.cmd)


def submit(ctx, cmd):
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(ctx.obj['OUTPUT'], 'logs/runner_{}'.format(r))
    workflow_runner = Runner(cmd=cmd, output_dp=ctx.obj['OUTPUT'], max_ppn=ctx.obj['PPN'], max_mem=ctx.obj['MEM'])
    workflow_runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])


@click.group()
@click.option('--output', '-o', default='metago_output')
@click.option('--flux/--no-flux', default=False)
@click.option('--account', '-a')
@click.option('--ppn', '-p', default=4)
@click.option('--mem', '-m', default='20000') # current limitation, only handles mb
@click.option('--walltime', '-w', default='2:00:00')
@click.pass_context
def cli(ctx, output, flux, account, ppn, mem, walltime):
    if not os.path.exists(output): os.makedirs(output)
    fp = os.path.dirname(os.path.abspath(__file__))
    environment = os.path.join(fp, 'dependencies/miniconda/bin/activate')
    qc_script = os.path.join(fp, 'workflow/quality_control.py')
     
    if flux:
        if not account: sys.exit('To attempt a submission to the flux cluster you need to supply an --account/-a')
        cmd = ' '.join(sys.argv)
        full_dp = os.path.dirname(os.path.abspath(__file__))
        runner_fp = os.path.join(full_dp, 'metaGO.py')
        qsub = 'qsub -N metaGO -A {} -q fluxm -l nodes=1:ppn={},mem={}mb,walltime={}'.format(
                                                                   account, ppn, mem, walltime)
        call('echo "source {} && python {}" | {}'.format(environment, cmd, qsub), shell=True)
        sys.exit('Launched command via Flux')
    ctx.obj['ENVIRONMENT'] = environment
    ctx.obj['QCSCRIPT'] = qc_script
    ctx.obj['OUTPUT'] = output
    ctx.obj['FLUX'] = flux
    ctx.obj['ACCOUNT'] = account
    ctx.obj['PPN'] = ppn
    ctx.obj['MEM'] = mem
    ctx.obj['WALLTIME'] = walltime


@cli.command()
@click.argument('sample_dp')
@click.pass_context
def sample_qc(ctx, sample_dp):
    """ Sample quality control against Illumina sample """
    click.echo('Sample QC called')
    cmd = 'source {} && python {} -o {} -p {} -m {} sample_qc {}'.format(ctx.obj['ENVIRONMENT'],
                                                                      ctx.obj['QCSCRIPT'],
                                                                      ctx.obj['OUTPUT'],
                                                                      ctx.obj['PPN'],
                                                                      ctx.obj['MEM'],
                                                                      sample_dp)
    submit(ctx=ctx, cmd=cmd)

@cli.command()
@click.argument('run_dp')
@click.pass_context
def run_qc(ctx, run_dp):
    """ Run quality control against entire Illumina run """
    click.echo('Run QC called')
    cmd = 'source {} && python {} -o {} -p {} -m {} run_qc {}'.format(ctx.obj['ENVIRONMENT'],
                                                                      ctx.obj['QCSCRIPT'],
                                                                      ctx.obj['OUTPUT'],
                                                                      ctx.obj['PPN'],
                                                                      ctx.obj['MEM'],
                                                                      run_dp)
    submit(ctx=ctx, cmd=cmd)


if __name__ == "__main__":
    cli(obj={})
