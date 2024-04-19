'''
This script divides the initial multifastas into individual fastas and gives the right description format to every seq
'''
import os
import shutil
from Bio import SeqIO

def genomes(multifasta):
    if os.path.isdir('genomes'):
        shutil.rmtree('genomes')
    os.mkdir('genomes')

    for seq in SeqIO.parse(multifasta, 'fasta'):
        description = seq.description[seq.description.find(' ')+1:]
        with open(f'genomes/{seq.id}.fasta', 'w', encoding='utf-8') as output:
            output.write(f'>{seq.id} {description}\n{seq.seq}\n')

def get_seqs(multifasta):
    '''
    OUT: dict[genome_id] = [biopythonseq1, biopythonseq2,..]
    '''
    seqs = {}
    for seq in SeqIO.parse(multifasta, 'fasta'):
        genome_id = seq.id[seq.id.find('|')+1:seq.id.find('.', seq.id.find('|'))+2]
        if genome_id not in seqs:
            seqs[genome_id] = []
        seqs[genome_id].append(seq)
    return seqs

def prot_orf(multifasta, dtype):
    folder = 'proteomes' if dtype == '_prot_' else 'orfeomes'

    if os.path.isdir(folder):
        shutil.rmtree(folder)
    os.mkdir(folder)

    seqs = get_seqs(multifasta) #dict[genome_id] = [biopythonseq1, biopythonseq2,..]
    for genome_id, seq_list in seqs.items():
        with open(f'{folder}/{genome_id}.fasta', 'w', encoding='utf-8') as output:
            for seq in seq_list:
                prot_id = seq.id[seq.id.find(dtype)+len(dtype):seq.id.find('.', seq.id.find(dtype))+2]
                description = seq.description[seq.description.find(' ')+1:]
                output.write(f'>{prot_id} {genome_id} {description}\n{seq.seq}\n')

def main():
    genomes('genomes.fasta')
    prot_orf('proteomes.fasta', '_prot_')
    prot_orf('orfeomes.fasta', '_cds_')