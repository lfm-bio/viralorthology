#!/usr/bin/python3
from modules.misc import combine_fastas
from tqdm import tqdm
import os

def get_ids():
    ids_file = open('ids.txt')
    ids = []
    for line in ids_file:
        ids.append(line.strip())
    ids_file.close()
    return ids

def download_fastas(ids):
    print('Downloading files...')
    for n, seqid in enumerate(tqdm(ids)):
        while True:
            os.system(f'efetch -db nuccore -id "{seqid}" -format fasta > new_genome{n}.fasta')
            os.system(f'efetch -db nuccore -id "{seqid}" -format fasta_cds_aa > new_proteome{n}.fasta') #proteome
            os.system(f'efetch -db nuccore -id "{seqid}" -format fasta_cds_na > new_orfeome{n}.fasta') #orfeome
            if os.path.exists(f'new_proteome{n}.fasta') and os.path.exists(f'new_orfeome{n}.fasta') and os.path.exists(f'new_genome{n}.fasta'):
                break

def make_mulfitastas():
    genomes = sorted([fasta for fasta in os.listdir(os.curdir) if fasta.startswith('new_genome')])
    proteomes = sorted([fasta for fasta in os.listdir(os.curdir) if fasta.startswith('new_proteome')])
    orfeomes = sorted([fasta for fasta in os.listdir(os.curdir) if fasta.startswith('new_orfeome')])
    combine_fastas(genomes, 'genomes.fasta', 'fasta')
    combine_fastas(proteomes, 'proteomes.fasta', 'fasta')
    combine_fastas(orfeomes, 'orfeomes.fasta', 'fasta')

def main(params = ''):
    ids = get_ids()
    download_fastas(ids)
    make_mulfitastas()