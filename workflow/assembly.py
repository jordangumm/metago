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
        
        To consider: estimating diversity (genome complexity) and depth to inform assembly parameter settings
        """
        pairs, is_paired = self.get_fastq_pairs(self.fastqs) # put in workflow instead of __init__ so print statements log
        fp = os.path.dirname(os.path.abspath(__file__))
        conda = os.path.join(fp, '../dependencies/miniconda/bin/activate')

        sample_output_dp = os.path.join(self.output_dp, self.sid)
        if not os.path.exists(sample_output_dp): os.makedirs(sample_output_dp)
        mem_available = self.getMemMb()
        if mem_available == 'unlimited':
            mem_available = self.max_mem
        else:
            mem_available = int(mem_available)



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

    workflow_runner = SampleAssembly(run_dp=run_dp, output_dp=output_dp, max_ppn=ppn, max_mem=mem)
    workflow_runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ppn, memMb=mem)


if __name__ == "__main__":
    runner()
