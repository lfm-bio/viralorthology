#!/usr/bin/python3
import os
import time
from tqdm import tqdm
from modules.misc import combine_fastas

def get_ids():
    ids_file = open('ids.txt', encoding='utf-8')
    ids = []
    for line in ids_file:
        ids.append(line.strip())
    ids_file.close()
    return ids

def download_fastas(ids):
    print('Downloading files...')
    for seqid in tqdm(ids):
        while True:
            err1 = os.system(f'efetch -db nuccore -id "{seqid}" -format fasta > {seqid}.genome')
            time.sleep(2)
            err2 = os.system(f'efetch -db nuccore -id "{seqid}" -format fasta_cds_aa > {seqid}.proteome')
            time.sleep(2)
            err3 = os.system(f'efetch -db nuccore -id "{seqid}" -format fasta_cds_na > {seqid}.orfeome')
            time.sleep(2)
            if sum([err1, err2, err3]) == 0:
                break

def make_mulfitastas():
    genomes = sorted([fasta for fasta in os.listdir(os.curdir) if fasta.endswith('.genome')])
    proteomes = sorted([fasta for fasta in os.listdir(os.curdir) if fasta.endswith('.proteome')])
    orfeomes = sorted([fasta for fasta in os.listdir(os.curdir) if fasta.endswith('.orfeome')])
    combine_fastas(genomes, 'genomes.fasta', 'fasta')
    combine_fastas(proteomes, 'proteomes.fasta', 'fasta')
    combine_fastas(orfeomes, 'orfeomes.fasta', 'fasta')

def main():
    ids = get_ids()
    download_fastas(ids)
    make_mulfitastas()