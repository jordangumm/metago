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

        self.max_ppn = max_ppn
        self.max_mem = max_mem

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
        return pairs, True

    def workflow(self):
        """ Quality Control Workflow
        
        To consider: integrating memory estimation to dynamically set requirements
        """
        pairs, is_paired = self.get_fastq_pairs(self.fastqs) # put in workflow instead of __init__ so print statements log
        fp = os.path.dirname(os.path.abspath(__file__))
        conda = os.path.join(fp, '../dependencies/miniconda/bin/activate')
        adapters = os.path.join(fp, '../dependencies/adapters.fa')
        sample_output_dp = os.path.join(self.output_dp, self.sid)
        if not os.path.exists(sample_output_dp): os.makedirs(sample_output_dp)
        mem_available = self.getMemMb()
        if mem_available == 'unlimited':
            mem_available = self.max_mem
        else:
            mem_available = int(mem_available)

        scheduled_tasks = []
        if is_paired:
            for pair in pairs.keys():
                pair_tasks = []

                pair_size = (os.path.getsize(pairs[pair]['r1']) >> 20) * 2 # in MB, both pairs
                print '\n\n\n{}\n\n\n'.format(pair_size)
                if pair_size > mem_available:
                    sys.exit('{}MB is not enough memory to interleave {}MB pair'.format(mem_available, pair_size))

                # Step 1. Interleave files to set up for next steps
                interleaved_fp = os.path.join(sample_output_dp, 'interleaved', '{}.fastq'.format(pair))
                if not os.path.exists(interleaved_fp):
                    cmd = 'source {} && reformat.sh t=4'.format(conda)
                    cmd += ' in1={} in2={} out={}'.format(pairs[pair]['r1'], pairs[pair]['r2'], interleaved_fp)
                    self.addTask("interleave_{}".format(pair), nCores=4, memMb=pair_size, command=cmd)
                    pair_tasks.append("interleave_{}".format(pair))

                # Step 2. Quality control by trimming adapters and low quality bases for mapping
                trimmed_fp = os.path.join(sample_output_dp, 'quality_controlled', '{}.fastq'.format(pair))
                if not os.path.exists(trimmed_fp):
                    cmd = 'source {} && bbduk.sh in={} out={} ref={}'.format(conda, interleaved_fp, trimmed_fp, adapters)
                    cmd += ' ktrim=r k=23 mink=11 hdist=1 tpe tbo t=1 qtrim=rl trimq=20 maq=20 interleaved=t'
                    self.addTask("trim_{}".format(pair), nCores=1, memMb=pair_size, command=cmd, dependencies=pair_tasks) 
                    pair_tasks.append("trim_{}".format(pair))

                # Step 3. Normalize to even coverage for assembly
                #normalized_fp = os.path.join(sample_output_dp, 'normalized', '{}.fastq'.format(pair))
                #if not os.path.exists(normalized_fp):
                #    cmd = 'source {} && bbnorm.sh t=4'.format(conda)
                #    cmd += ' in={} out={} target=100 min=5'.format(trimmed_fp, normalized_fp)
                #    self.addTask("normalize_{}".format(pair), nCores=4, memMb=pair_size*4, command=cmd, dependencies=pair_tasks)
                #    pair_tasks.append("normalize_{}".format(pair))
                
                scheduled_tasks += pair_tasks
        else:
            sys.exit('Not Implemented: single paired QC')

        sid = pair.split('_')[0]
        merged_fp = os.path.join(sample_output_dp, '{}.fastq'.format(sid))
        if not os.path.exists(merged_fp):
            cmd = 'source {} && cat {}/* > {}'.format(conda, os.path.dirname(os.path.abspath(trimmed_fp)), merged_fp)
            self.addTask("join_{}".format(pair), nCores=1, memMb=2000, command=cmd, dependencies=scheduled_tasks)
            scheduled_tasks.append("join_{}".format(pair))

        normalized_fp = os.path.join(sample_output_dp, 'normalized', '{}.fastq'.format(sid))
        if not os.path.exists(normalized_fp):
            cmd = 'source {} && bbnorm.sh t=4'.format(conda)
            cmd += ' in={} out={} target=100 min=5'.format(merged_fp, normalized_fp)
            fastq_size = (os.path.getsize(merged_fp) >> 20) * 4 # in MB adjusted for bbnorm requirements
            num_cores = self.max_mem / fastq_size # take up cores based on ratio of memory to consume
            self.addTask("normalize_{}".format(pair), nCores=num_cores, memMb=fastq_size, command=cmd, dependencies=pair_tasks)
            scheduled_tasks.append("normalize_{}".format(pair))


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
                if 'fastq' not in fastq: continue
                fastq_fp = os.path.join(sample_dp, fastq)
                sample_fastqs.append(fastq_fp)
            sample_qc_runner = SampleQualityControl(sid=sample, fastqs=sample_fastqs, output_dp=self.output_dp,
                                                                     max_ppn=self.max_ppn, max_mem=self.max_mem)
            self.addWorkflowTask(label=sample, workflowRunnerInstance=sample_qc_runner)
            break


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
    log_output_dp = os.path.join(output_dp, 'logs', 'quality_control')

    workflow_runner = RunQualityControl(run_dp=run_dp, output_dp=output_dp, max_ppn=ppn, max_mem=mem)
    workflow_runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ppn, memMb=mem)


if __name__ == "__main__":
    runner()
