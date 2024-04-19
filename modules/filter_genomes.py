'''
This script filters similar genomes, and moves them to another folder (/filtered)
'''
import os
from modules.misc import get_file_list
from modules.misc import delete_tmp_files
from modules import commands

def makeblastdb():
    fastas = get_file_list()
    for fasta in fastas:
        commands.make_blastdb_nt(fasta)

def are_similar():
    similar = False
    with open('test', encoding='utf-8') as result:
        for line in result:
            if '# 0 hits found' in line:
                break
            if '@' in line:
                q_cov =  int(line.strip().split('@')[0])
                ident = float(line.strip().split('@')[1])
                if q_cov > 90 and ident > 95: # PARAMS TO DETERMINE SIMILITUDE
                    similar = True
                    break
    os.remove('test')
    return similar

def filter_genomes():
    fastas = get_file_list()
    print('\nFiltering genomes')
    for n, f1 in enumerate(fastas):
        porc = round((n+1) / len(fastas) * 100, 2)
        print(f'\r{porc}%', end='')
        for f2 in fastas:
            if f2 == f1:
                continue
            commands.blast_filter_genomes(f1, f2)
            if are_similar():
                if f2 in fastas:
                    fastas.pop(fastas.index(f2))
    print()
    return fastas

def move_files(fastas):
    os.chdir('..')
    paths = ['genomes', 'proteomes', 'orfeomes']
    for path in paths:
        os.chdir(path)
        os.mkdir('filtered')
        fastas_all = get_file_list()
        for fasta in fastas_all:
            if fasta not in fastas:
                os.replace(fasta, f'filtered/{fasta}')
        os.chdir('..')

def main():
    os.chdir('genomes')

    makeblastdb()
    fastas = filter_genomes()
    delete_tmp_files(['.nhr', '.nin', '.nsq'])
    move_files(fastas) #final chdir in func

if __name__ == '__main__':
    main()