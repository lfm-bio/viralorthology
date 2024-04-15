import os
from Bio import SeqIO
from modules.misc import get_file_list
from modules.misc import get_ordered_files

def makeDB(file_ext, db_name):
    proteomes = get_file_list(file_ext)
    proteomes = get_ordered_files(proteomes)
    with open(f'../{db_name}', 'w', encoding='utf-8') as output:
        for proteome in proteomes:
            for seq in SeqIO.parse(proteome, 'fasta'):
                output.write(seq.format('fasta-2line'))
    output.close()

def main():
    os.chdir('proteomes')

    makeDB('.fasta', 'protDB.db')
    makeDB('.OFprot', 'protDB_OF.db')

    os.chdir('..')