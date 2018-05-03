import sys
import os
import random
import string

import click
from pyflux import FluxWorkflowRunner
from subprocess import call

class SampleQualityControl(FluxWorkflowRunner):
    """ Follows JGI Data Preprocessing Guide

    https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/data-preprocessing/
    """
    def __init__(self, sid, fastqs, output_dp, max_ppn, max_mem):
        self.sid = sid
        self.fastqs = fastqs
        self.output_dp = output_dp

        self.max_ppn = int(max_ppn)
        self.max_mem = int(max_mem)

    def get_fastq_pairs(self, fastqs):
        """ Return lists of fastq pairs

        Uses Illumina sample naming conventions to identify pairs.
        Returns original list if no pairs found.
        """
        pairs = {}
        for fastq in fastqs:
            if '_R1_' in fastq:
                pair_name = '_'.join(fastq.split('/')[-1].split('_R1_')).replace('.fastq','').replace('.gz','')
                pairs[pair_name] = {}
                if fastq.replace('_R1_', '_R2_') in fastqs:
                    pairs[pair_name]['r1'] = fastq
                    pairs[pair_name]['r2'] = fastq.replace('_R1_', '_R2_')
                else: # not paired
                    print("NOT PAIRED: {}".format(fastq))
                    return fastqs, False
            if '_1.fa' in fastq:
                pair_name = '_'.join(fastq.split('/')[-1].split('_1.fa')).replace('stq','').replace('.gz','').replace('_','')
                pairs[pair_name] = {}
                if fastq.replace('_1.fa', '_2.fa') in fastqs:
                    pairs[pair_name]['r1'] = fastq
                    pairs[pair_name]['r2'] = fastq.replace('_1.fa', '_2.fa')
                else: # not paired
                    print("NOT PAIRED: {}".format(fastq))
                    return fastqs, False
        return pairs, True

    def workflow(self):
        """ Quality Control Workflow
        
        To consider: integrating memory estimation to dynamically set requirements
        """
        pairs, is_paired = self.get_fastq_pairs(self.fastqs) # put in workflow instead of __init__ so print statements log
        if len(pairs) == 0: sys.exit("No pairs found")
        fp = os.path.dirname(os.path.abspath(__file__))
        conda = os.path.join(fp, '../dependencies/miniconda/bin/activate')
        adapters = os.path.join(fp, '../dependencies/adapters.fa')
        sample_output_dp = os.path.join(self.output_dp, self.sid)
        if not os.path.exists(sample_output_dp): os.makedirs(sample_output_dp)

        scheduled_tasks = []
        if is_paired:
            for pair in pairs.keys():
                pair_tasks = []

                pair_size = (os.path.getsize(pairs[pair]['r1']) >> 20) * 2 # in MB, both pairs
                if pair_size > self.max_mem:
                    sys.exit('{}MB is not enough memory to interleave {}MB pair'.format(self.max_mem, pair_size))
                if pair_size == 0: continue # bad sample

                # Step 1. Interleave files to set up for next steps
                interleaved_fp = os.path.join(sample_output_dp, 'interleaved', '{}.fastq'.format(pair))
                if not os.path.exists(interleaved_fp):
                    cmd = 'source {} && reformat.sh t=4'.format(conda)
                    cmd += ' in1={} in2={} out={}'.format(pairs[pair]['r1'], pairs[pair]['r2'], interleaved_fp)
                    self.addTask("interleave_{}".format(pair), nCores=4, memMb=pair_size*2, command=cmd)
                    pair_tasks.append("interleave_{}".format(pair))

                # Step 2. Quality control by trimming adapters and low quality bases for mapping
                trimmed_fp = os.path.join(sample_output_dp, 'quality_controlled', '{}.fastq'.format(pair))
                if not os.path.exists(trimmed_fp):
                    cmd = 'source {} && bbduk.sh in={} out={} ref={}'.format(conda, interleaved_fp, trimmed_fp, adapters)
                    cmd += ' ktrim=r k=23 mink=11 hdist=1 tpe tbo t=4 qtrim=rl trimq=20 maq=20 interleaved=t'
                    self.addTask("trim_{}".format(pair), nCores=4, memMb=pair_size*2, command=cmd, dependencies=pair_tasks) 
                    pair_tasks.append("trim_{}".format(pair))
                scheduled_tasks += pair_tasks

        else:
            sys.exit('Not Implemented: single end QC')

        merged_fp = os.path.join(sample_output_dp, '{}.fastq'.format(self.sid))
        if not os.path.exists(merged_fp):
            cmd = 'source {} && cat {}/* > {}'.format(conda, os.path.dirname(os.path.abspath(trimmed_fp)), merged_fp)
            self.addTask("join_{}".format(pair), nCores=1, memMb=2000, command=cmd, dependencies=scheduled_tasks)
            scheduled_tasks.append("join_{}".format(pair))


class RunQualityControl(FluxWorkflowRunner):
    def __init__(self, run_dp, output_dp, max_ppn, max_mem):
        self.run_dp = run_dp
        self.output_dp = output_dp

        self.max_ppn = max_ppn
        self.max_mem = max_mem

    def workflow(self):
        for sample in os.listdir(self.run_dp):
            sample_dp = os.path.join(self.run_dp, sample)
            if not os.path.isdir(sample_dp) or not 'Sample_' in sample: continue
            sample_fastqs = []
            for fastq in os.listdir(sample_dp):
                if '.fastq' not in fastq and not '.fa' in fastq: continue
                fastq_fp = os.path.join(sample_dp, fastq)
                sample_fastqs.append(fastq_fp)

            sample_qc_runner = SampleQualityControl(sid=sample, fastqs=sample_fastqs, output_dp=self.output_dp,
                                                                     max_ppn=self.max_ppn, max_mem=self.max_mem)
            self.addWorkflowTask(label=sample, workflowRunnerInstance=sample_qc_runner)


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
@click.pass_context
def run_qc(ctx, run_dp):
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output = os.path.join(ctx.obj['OUTPUT'], 'logs/run_qc_{}'.format(r))
    runner = RunQualityControl(run_dp=run_dp,
                               output_dp=ctx.obj['OUTPUT'],
                               max_ppn=ctx.obj['PPN'],
                               max_mem=ctx.obj['MEM'])
    return_code = runner.run(mode='local', dataDirRoot=log_output, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])
    if return_code != 0: sys.exit('Non-Zero Exit Code in run_qc')
    return return_code


@cli.command()
@click.argument('sample_dp')
@click.pass_context
def sample_qc(ctx, sample_dp):
    sid = sample_dp.split('/')[-1]
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output = os.path.join(ctx.obj['OUTPUT'], sid, 'qc_{}'.format(r))

    sample_fastqs = []
    for fastq in os.listdir(sample_dp):
        if '.fastq' not in fastq and '.fa' not in fastq: continue
        fastq_fp = os.path.join(sample_dp, fastq)
        sample_fastqs.append(fastq_fp) 
    if len(sample_fastqs) == 0: sys.exit('No fastq files in {}'.format(sample_dp))

    runner = SampleQualityControl(sid=sid,
                                  fastqs=sample_fastqs,
                                  output_dp=ctx.obj['OUTPUT'],
                                  max_ppn=ctx.obj['PPN'],
                                  max_mem=ctx.obj['MEM'])
    return_code = runner.run(mode='local', dataDirRoot=log_output, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])
    if return_code != 0: sys.exit('Non-Zero Exit Code in run_qc')
    return return_code


if __name__ == "__main__":
    cli(obj={})
