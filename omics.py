import click
import os

import screed
import pandas as pd


@click.group()
def cli():
    pass


@cli.command()
@click.option('--fasta_fp', '-f', help='path to fasta file to split', required=True)
@click.option('--ref_fp', '-r', help='path to reference fasta file', required=True)
@click.option('--output_dp', '-o', help='path to directory to split fasta to', required=True)
def split(fasta_fp, ref_fp, output_dp):
    if not os.path.exists(output_dp): os.makedirs(output_dp)
    fasta_db = screed.read_fasta_sequences(fasta_fp)

    n = 0
    first_entry = True
    for i, seq in enumerate(fasta_db):
        if i == 0:
            output = open(os.path.join(output_dp, 'tmp_{}.fasta'.format(n)), 'w+')
        elif i % split_num == 0:
            n += 1
            output.close()
            first_entry = True
            output = open(os.path.join(output_dp, 'tmp_{}.fasta'.format(n)), 'w+')
            
        if first_entry:
            output.write('>{}\n{}'.format(fasta_db[seq].name, str(fasta_db[seq].sequence)))
            first_entry = False
        else:
            output.write('\n>{}\n{}'.format(fasta_db[seq].name, str(fasta_db[seq].sequence)))
    output.close()
    


if __name__ == "__main__":
    cli()
