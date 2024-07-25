import os
import shutil
from Bio import SeqIO

def split_genomes(fasta):
    if os.path.isdir('genomes'):
        shutil.rmtree('genomes')
    os.mkdir('genomes')

    for seq in SeqIO.parse(fasta, 'fasta'):
        description = seq.description[seq.description.find(' ')+1:]
        seq_id = seq.id[:seq.id.find('.')]
        with open(f'genomes/{seq_id}.fasta', 'w', encoding='utf-8') as output:
            output.write(f'>{seq_id} {description}\n{seq.seq}\n')

def split_prots_orfs(folder):

    def get_seqs(fasta):
        '''
        OUT: dict[genome_id] = [biopythonseq1, biopythonseq2,..]
        '''
        seqs = {}
        for seq in SeqIO.parse(fasta, 'fasta'):
            genome_id = seq.id[seq.id.find('|')+1:seq.id.find('.', seq.id.find('|'))]
            if genome_id not in seqs:
                seqs[genome_id] = []
            seqs[genome_id].append(seq)
        return seqs

    dtype = '_prot_' if folder == 'proteomes' else '_cds_'
    fasta = folder + '.fasta'

    if os.path.isdir(folder):
        shutil.rmtree(folder)
    os.mkdir(folder)

    seqs = get_seqs(fasta) #dict[genome_id] = [biopythonseq1, biopythonseq2,..]
    for genome_id, seq_list in seqs.items():
        with open(f'{folder}/{genome_id}.fasta', 'w', encoding='utf-8') as output:
            for seq in seq_list:
                prot_id = seq.id[seq.id.find(dtype)+len(dtype):seq.id.find('.', seq.id.find(dtype))]
                description = seq.description[seq.description.find(' ')+1:]
                output.write(f'>{prot_id} {genome_id} {description}\n{seq.seq}\n')

def main():
    split_genomes('genomes.fasta')
    split_prots_orfs('proteomes')
    split_prots_orfs('orfeomes')
