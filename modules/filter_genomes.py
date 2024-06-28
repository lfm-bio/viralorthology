import os
from tqdm import tqdm
from modules.misc import get_file_list
from modules.misc import delete_tmp_files
from modules import commands

def filter_genomes():

    def genomes_are_similar():
        similar = False
        with open('blast.results', encoding='utf-8') as result:
            for line in result:
                if '# 0 hits found' in line:
                    break
                if '@' in line:
                    q_cov =  int(line.strip().split('@')[0])
                    ident = float(line.strip().split('@')[1])
                    if q_cov > 90 and ident > 95: # PARAMS TO DETERMINE SIMILITUDE
                        similar = True
                        break
        os.remove('blast.results')
        return similar

    fastas = get_file_list()
    for fasta in fastas:
        commands.make_blastdb_nt(fasta)

    print('Filtering genomes...')
    for f1 in tqdm(fastas):
        for f2 in fastas:
            if f2 == f1:
                continue
            commands.blast_filter_genomes(f1, f2)
            if genomes_are_similar() and f2 in fastas:
                fastas.pop(fastas.index(f2))

    return fastas

def move_files(fastas):
    os.chdir('..')
    folders = ['genomes', 'proteomes', 'orfeomes']
    for folder in folders:
        os.chdir(folder)
        fastas_all = get_file_list()
        if len(fastas_all) != len(fastas):
            os.mkdir('filtered')
            for fasta in fastas_all:
                if fasta not in fastas:
                    os.replace(fasta, f'filtered/{fasta}')
        os.chdir('..')

def main():
    os.chdir('genomes')

    fastas = filter_genomes()
    delete_tmp_files(['.nhr', '.nin', '.nsq'])
    move_files(fastas)

    #final chdir in move_files()
