import sys
import os
import random
import string

import click
from pyflux import FluxWorkflowRunner
from subprocess import call

import logging

class SampleQualityControl(FluxWorkflowRunner):
    """ Follows JGI Data Preprocessing Guide

    https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/data-preprocessing/
    """
    def __init__(self, sid, fastqs, output_dp, max_ppn, max_mem, overwrite):
        self.sid = sid
        self.fastqs = fastqs
        self.output_dp = output_dp

        self.max_ppn = int(max_ppn)
        self.max_mem = int(max_mem)
        self.overwrite = overwrite

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
        fastqs, is_paired = self.get_fastq_pairs(self.fastqs) # put in workflow instead of __init__ so print statements log
        print fastqs
        if len(fastqs) == 0:
            print("No fastq found: is the post-QC output?")
            return
        fp = os.path.dirname(os.path.abspath(__file__))
        conda = os.path.join(fp, '../dependencies/miniconda/bin/activate')
        adapters = os.path.join(fp, '../dependencies/adapters.fa')
        sample_output_dp = os.path.join(self.output_dp, self.sid)
        if not os.path.exists(sample_output_dp):
            os.makedirs(sample_output_dp)
        else:
            logging.warning('[ERROR]: bbduk.sh requires output to be fresh: force using --overwrite option')
            sys.exit('bbduk.sh requires output to be fresh: force using --overwrite option')

        scheduled_tasks = []
        if is_paired:
            for pair in fastqs.keys():
                pair_tasks = []

                r1_size = (os.path.getsize(fastqs[pair]['r1']) >> 20) # in MB, forward strands
                r2_size = (os.path.getsize(fastqs[pair]['r2']) >> 20) # in MB, reverse strands
                pair_size = r1_size + r2_size
                if pair_size > self.max_mem:
                    logging.warning('Possible failure ahead; low memory supplied for handling {}MB fastq pair'.format(pair_size))
                rsizediff = abs(r2_size-r1_size)
                if rsizediff > r2_size*0.10 or rsizediff > r1_size*0.10:
                    logging.warning('Possible issue with pair; {} greater than 10% different in size'.format(pair))
                if r1_size == 0:
                    logging.warning('Possible failed sequencing; {} is less than 1MB in size'.format(fastqs[pair]['r1']))
                if r2_size == 0:
                    logging.warning('Possible failed sequencing; {} is less than 1MB in size'.format(fastqs[pair]['r2']))
                if pair_size == 0: pair_size = 1

                # Step 0. Correct any issues? -- so far issues have been from not-completely-downloaded fastq files 
                # Step 1. Interleave files to set up for next steps
                interleaved_fp = os.path.join(sample_output_dp, 'interleaved', '{}.fastq'.format(pair))
                if not os.path.exists(interleaved_fp):
                    cmd = 'source {} && reformat.sh -Xmx{}m tossbrokenreads=t t=4'.format(conda, pair_size*2)
                    cmd += ' in1={} in2={} out={}'.format(fastqs[pair]['r1'], fastqs[pair]['r2'], interleaved_fp)
                    self.addTask("interleave_{}".format(pair), nCores=4, memMb=pair_size*2, command=cmd)
                    pair_tasks.append("interleave_{}".format(pair))

                # Step 2. Quality control by trimming adapters and low quality bases for mapping
                trimmed_fp = os.path.join(sample_output_dp, 'quality_controlled', '{}.fastq'.format(pair))
                stats_fp = os.path.join(sample_output_dp, 'quality_controlled', 'stats.txt')
                cmd = 'source {} && bbduk.sh -Xmx{}m in={} out={} ref={} stats={}'.format(conda, pair_size*2, interleaved_fp, trimmed_fp, adapters, stats_fp)
                cmd += ' ktrim=r k=23 mink=11 hdist=1 tpe tbo t=4 qtrim=rl trimq=20 maq=20 interleaved=t minlen=70'
                if self.overwrite: cmd += ' overwrite=t'
                self.addTask("trim_{}".format(pair), nCores=4, memMb=pair_size*2, command=cmd, dependencies=pair_tasks) 
                pair_tasks.append("trim_{}".format(pair))
                scheduled_tasks += pair_tasks

        else:
            for fastq in fastqs:
                fastq_tasks = []
                fastq_size = (os.path.getsize(fastq) >> 20) # in MB
                if fastq_size > self.max_mem:
                    logging.warning('Possible failure ahead; low memory supplied for handling {}MB fastq pair'.format(fastq_size))
                if fastq_size == 0:
                    logging.warning('Possible failed sequencing; {} is less than 1MB in size'.format(fastq))
                    fastq_size = 1
                # Step 1. Quality control by trimming adapters and low quality bases for mapping
                fastq_name = fastq.split('/')[-1].replace('.fastq','').replace('.gz','')
                trimmed_fp = os.path.join(sample_output_dp, 'quality_controlled', '{}.fastq'.format(fastq_name))
                cmd = 'source {} && bbduk.sh -Xmx{}m in={} out={} ref={} t=4'.format(conda, fastq_size*2, fastq, trimmed_fp, adapters)
                if self.overwrite: cmd += ' overwrite=t'
                #cmd += ' t=4 trimq=20 maq=20 minlen=70'
                #cmd += ' ktrim=r k=23 mink=11 hdist=1 t=4 qtrim=rl trimq=20 maq=20 minlen=70'
                self.addTask("trim_{}".format(fastq_name), nCores=4, memMb=fastq_size*2, command=cmd)
                scheduled_tasks.append("trim_{}".format(fastq_name))

        merged_fp = os.path.join(sample_output_dp, '{}.fastq'.format(self.sid))
        if not os.path.exists(merged_fp):
            cmd = 'source {} && cat {}/* > {}'.format(conda, os.path.join(sample_output_dp, 'quality_controlled'), merged_fp)
            self.addTask("join_fastq", nCores=1, memMb=2000, command=cmd, dependencies=scheduled_tasks)
            scheduled_tasks.append("join_fastq")


class RunQualityControl(FluxWorkflowRunner):
    def __init__(self, run_dp, output_dp, max_ppn, max_mem, overwrite):
        self.run_dp = run_dp
        self.output_dp = output_dp

        self.max_ppn = max_ppn
        self.max_mem = max_mem
        self.overwrite = overwrite

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
                                           max_ppn=self.max_ppn, max_mem=self.max_mem, overwrite=self.overwrite)
            self.addWorkflowTask(label=sample, workflowRunnerInstance=sample_qc_runner)


def run_qc(ctx, run_dp, overwrite):
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output = os.path.join(ctx.obj['OUTPUT'], 'logs/run_qc_{}'.format(r))
    runner = RunQualityControl(run_dp=run_dp,
                               output_dp=ctx.obj['OUTPUT'],
                               max_ppn=ctx.obj['PPN'],
                               max_mem=ctx.obj['MEM'],
                               overwrite=overwrite)
    return_code = runner.run(mode='local', dataDirRoot=log_output, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])
    if return_code != 0: sys.exit('Non-Zero Exit Code in run_qc')
    return return_code


def sample_qc(ctx, sample_dp, overwrite):
    sid = sample_dp.rstrip('/').split('/')[-1]
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
                                  max_mem=ctx.obj['MEM'],
                                  overwrite=overwrite)
    return_code = runner.run(mode='local', dataDirRoot=log_output, nCores=ctx.obj['PPN'], memMb=ctx.obj['MEM'])
    if return_code != 0:
        logger.warning('[ERROR]: Non-Zero Exit Code in run_qc')
        sys.exit('Non-Zero Exit Code in run_qc')
    return return_code


@click.command()
@click.option('--run', help='directory path to Illumina run containing sample directories')
@click.option('--sample', help='directory path to Illumina sample containing fastq(s)')
@click.option('--fastq', help='path(s) to single sample single-end fastq(s)', multiple=True)
@click.option('--r1', help='comma seperated list of single sample forward strand fastq(s)')
@click.option('--r2', help='comma seperated list of single sample reverse strand fastq(s)')
@click.option('--overwrite/--not-overwrite', help='force bbduk.sh to overwrite files if needed', default=False)
@click.pass_context
def qc(ctx, run, sample, fastq, r1, r2, overwrite):
    """ Quality control fastq file(s) """
    empty_output = True
    if os.path.exists(ctx.obj['OUTPUT']):
        if len(os.listdir(ctx.obj['OUTPUT'])) > 0:
            empty_output = False
    logging.basicConfig(filename=os.path.join(ctx.obj['OUTPUT'], 'error.log'))
    
    if not empty_output and not overwrite:
        logging.warning('[ERROR]: bbduk.sh requires output to be fresh: force by using the --overwrite option')
        sys.exit('bbduk.sh requires output to be fresh: force by using the --overwrite option')

    if not run and not sample and not fastq and not r1 and not r2:
        logging.warning('[ERROR]: No option supplied to specify fastq(s) for quality control!')
        sys.exit('No option supplied to specify fastq(s) for quality control!')

    if run: run_qc(ctx, run, overwrite)
    elif sample: sample_qc(ctx, sample, overwrite)
    elif fastq:
        logging.warning('[ERROR]: Single end manually defined fastq method not implemented')
        sys.exit('Single end manually defined fastq method not implemented')
    elif r1 or r2:
        logging.warning('[ERROR]: Paired end manually defined fastq method not implemented')
        sys.exit('Paired end manually defined fastq method not implemented')
    else:
        logging.warning('[ERROR]: Unexpected failure: debug the logs!')
        sys.exit('Unexpected failure: debug')
