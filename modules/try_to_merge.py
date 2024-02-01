import os
import sys
from modules.misc import combine_fastas

def check_compatibility(fasta_list):
    genomes = []
    for fasta in fasta_list:
        fasta = open(fasta)
        for line in fasta:
            if line.startswith('>'):
                genome = line.strip().split()[1]
                genomes.append(genome)
        fasta.close()
    if len(genomes) == len(set(genomes)):
        return True
    else:
        print('Incompatible groups')
        print('Repetead genomes: ')
        for genome in set(genomes):
            if genomes.count(genome) == 2:
                print(genome)

        return False

def try_to_merge(fasta_list):
    if check_compatibility(fasta_list):
        combine_fastas(fasta_list, fasta_list[0])
        print(f'Groups merged. New name: {fasta_list[0]}')

def main(fasta_list):
    os.chdir('orthology_groups')

    try_to_merge(fasta_list)

    os.chdir('..')