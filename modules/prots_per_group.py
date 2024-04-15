import os
from modules.misc import get_file_list
from modules.misc import get_genomes_fasta

def get_groups():
    '''
    OUT: dict[group_name] = [genomes]
    '''
    f_groups = open('groups.csv', encoding='utf-8')
    groups = {}
    for line in f_groups:
        if line.startswith('>'):
            group_name = line.strip()[1:]
            continue
        line = sorted(line.strip().split(','))
        groups[group_name] = line
    f_groups.close()
    return groups

def get_genomes_per_fasta():
    fastas = get_file_list()

    genomes_per_fasta = {}
    for fasta in fastas:
        genomes_per_fasta[fasta] = sorted(get_genomes_fasta(fasta))

    return genomes_per_fasta


def submain(groups):
    output = open('../groups_per_orthologyGroup.txt', 'w', encoding='utf-8')
    genomes_per_fasta = get_genomes_per_fasta()

    for fasta in genomes_per_fasta:
        fasta_groups = []
        n_full_groups = 0
        for group in groups:
            if set(groups[group]).issubset(set(genomes_per_fasta[fasta])): #every genome of that group has a protein in that fasta
                fasta_groups.append(group.upper())
                n_full_groups += 1
            else:
                for genome in genomes_per_fasta[fasta]:
                    if genome in groups[group]:
                        fasta_groups.append(group.lower())
                        break

        if fasta_groups:
            if n_full_groups == len(groups): # coregene
                output.write(f'{fasta.replace(".fasta", ""):40s} CoreGene\n')
            else:
                fasta_groups.sort()
                output.write(f'{fasta.replace(".fasta", ""):40s} {(" ").join((fasta_groups))}\n')

    output.close()

def main():

    groups = get_groups()

    os.chdir('orthology_groups')

    submain(groups)

    os.chdir('..')