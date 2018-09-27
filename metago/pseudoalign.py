import sys
import os
import random
import string

import click
from pyflux import FluxWorkflowRunner
from subprocess import call


class SamplePseudoalignment(FluxWorkflowRunner):
    """ Coordinates sample read assignment counting """
    def __init__(self, fastq, reference, output, max_ppn, max_mem, paired):
        self.fastq = fastq
        self.reference = reference
        self.output = output
        self.max_ppn = max_ppn
        self.max_mem = max_mem
        self.paired = paired

    def workflow(self):
        fp = os.path.dirname(os.path.abspath(__file__))
        env = os.path.join(fp, '../dependencies/miniconda/bin/activate')
        output_dp = os.path.join(self.output, os.path.basename(self.fastq).split('.')[0])
        if not os.path.exists(output_dp): os.makedirs(output_dp)

        submitted_cmds = []

        # 1. generate kallisto kmer index for fast quantification
        kallisto_index = '{}.index'.format(self.reference)
        if not os.path.exists(kallisto_index):
            cmd = 'source {} && '.format(env)
            cmd += 'kallisto index -i {} {}'.format(kallisto_index, self.reference)
            self.addTask('kallisto_index', nCores=self.max_ppn, memMb=self.max_mem, command=cmd, dependencies=submitted_cmds)
            submitted_cmds.append('kallisto_index')

        if self.paired:
            # 2. deinterleave fastq so kallisto can quantify pseudoalignments
            deinterleaved_dp = os.path.join(os.path.dirname(os.path.abspath(self.fastq)), 'deinterleaved')
            deinterleaved_r1_fp = os.path.join(deinterleaved_dp, os.path.basename(self.fastq).replace('.fastq', '_R1_final.fastq'))
            deinterleaved_r2_fp = os.path.join(deinterleaved_dp, os.path.basename(self.fastq).replace('.fastq', '_R2_final.fastq'))
            if not os.path.exists(deinterleaved_r1_fp) or not os.path.exists(deinterleaved_r2_fp):
                cmd = 'source {} && '.format(env)
                cmd += 'reformat.sh in={} out1={} out2={}'.format(self.fastq, deinterleaved_r1_fp, deinterleaved_r2_fp)
                self.addTask('deinterleave', nCores=self.max_ppn, memMb=self.max_mem, command=cmd, dependencies=submitted_cmds)
                submitted_cmds.append('deinterleave')

            # 3. quantify pseudoalignment counts
            sam_fp = os.path.join(output_dp, '{}.sam'.format(self.reference.split('/')[-1].split('.')[0]))        
            if not os.path.exists(sam_fp):
                cmd = 'source {} && '.format(env) 
                cmd += 'kallisto quant --pseudobam -i {} -o {} -b 40 --bias --single-overhang --fr-stranded {} {} -t {} > {}'.format(
                                     kallisto_index, output_dp, deinterleaved_r1_fp, deinterleaved_r2_fp, self.max_ppn, sam_fp)
                self.addTask('kallisto', nCores=self.max_ppn, memMb=self.max_mem, command=cmd, dependencies=submitted_cmds)
                submitted_cmds.append('kallisto')
        else:
            # 2. quantify pseudoalignment counts
            sam_fp = os.path.join(output_dp, '{}.sam'.format(self.reference.split('/')[-1].split('.')[0]))
            if not os.path.exists(sam_fp):
                cmd = 'source {} && '.format(env)
                cmd += 'kallisto quant --pseudobam -i {} -o {} -b 40 --bias --single-overhang --single -l 150 {} -t {} > {}'.format(
                                     kallisto_index, output_dp, self.fastq, self.max_ppn, sam_fp)
                self.addTask('kallisto', nCores=self.max_ppn, memMb=self.max_mem, command=cmd, dependencies=submitted_cmds)
                submitted_cmds.append('kallisto')


class RunPseudoalignment(FluxWorkflowRunner):
    """ Coordinate multiple sample read assignment counting from a run directory
    
    Assumes samples have already been processed by quality control step into an interleaved state.
    """
    def __init__(self, run_dp, reference, output, max_ppn, max_mem, paired):
        self.run_dp = run_dp
        self.reference = reference
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
            sample_runner = SamplePseudoalignment(fastq=sample_fastq,
                                                  reference=self.reference,
                                                  output=self.output,
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
@click.argument('run_dp')
@click.option('--reference', '-r', required=True)
@click.option('--paired/--not-paired', default=False)
@click.pass_context
def run_pseudoalign(ctx, run_dp, reference, paired):
    """ Coordinate assignment of reads from run samples to reference """
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output = os.path.join(ctx.obj['OUTPUT'], 'logs/run_pseudoalign_{}'.format(r))
    runner = RunPseudoalignment(run_dp=run_dp,
                                reference=reference,
                                output=ctx.obj['OUTPUT'],
                                max_ppn=ctx.obj['PPN'],
                                max_mem=ctx.obj['MEM'],
                                paired=paired)
    return_code = runner.run(mode='local', dataDirRoot=log_output, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])
    if return_code != 0: sys.exit('Non-Zero Exit Code in run_pseudoalign')
    return return_code


@cli.command()
@click.argument('sample_fp')
@click.option('--reference', '-r', required=True)
@click.option('--paired/--not-paired', default=False)
@click.pass_context
def sample_pseudoalign(ctx, sample_fp, reference, paired):
    """ Assign reads from sample fastqs/fastas to reference """
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output = os.path.join(ctx.obj['OUTPUT'], 'logs/sample_pseudoalign_{}'.format(r))

    runner = SamplePseudoalignment(fastq=sample_fp,
                                   reference=reference,
                                   output=ctx.obj['OUTPUT'],
                                   max_ppn=ctx.obj['PPN'],
                                   max_mem=ctx.obj['MEM'],
                                   paired=paired)
    return_code = runner.run(mode='local', dataDirRoot=log_output, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])
    if return_code != 0: sys.exit('Non-Zero Exit Code in sample_pseudoalign')
    return return_code


if __name__ == "__main__":
    cli(obj={})
        
