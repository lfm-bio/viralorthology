import os
from modules.misc import get_file_list
from modules.misc import get_genomes_fasta
from modules.misc import combine_fastas

def check_compatibility(fastas):
    genomes = []
    for fasta in fastas:
        genomes += get_genomes_fasta(fasta)
    if len(genomes) == len(set(genomes)):
        return True
    else:
        return False

def ask_user(fastas):
    while True:
        print('The following files can be combined:')
        print((' / ').join(fastas))
        yes_no = input('Should I merge them? [y/n]: ')
        if yes_no in ['Y', 'y']:
            return True
        elif yes_no in ['N', 'n']:
            return False

def submain():
    fastas = get_file_list()
    skip = []
    for fasta in fastas:
        if 'hypothetical-protein' in fasta:
            continue
        group = [fasta]
        for f2 in fastas:
            if f2 == fasta or '-' not in f2 or f2 in skip:
                continue
            if fasta.replace('.fasta', '') in f2[:f2.rfind('-')]:
                group.append(f2)
        if len(group) > 1:
            if check_compatibility(group):
                if ask_user(group):
                    skip += group
                    combine_fastas(group, group[0])

def main():
    os.chdir('orthology_groups')

    submain()

    os.chdir('..')