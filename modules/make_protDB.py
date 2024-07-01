import os
from Bio import SeqIO
from modules.misc import get_file_list
from modules.misc import get_ordered_files

def make_db(file_ext, db_name):
    fastas = get_file_list(file_ext)
    fastas = get_ordered_files(fastas)
    with open(f'../{db_name}', 'w', encoding='utf-8') as output:
        for proteome in fastas:
            for seq in SeqIO.parse(proteome, 'fasta'):
                output.write(seq.format('fasta'))
    output.close()

def main():
    os.chdir('proteomes')

    make_db('.fasta', 'protDB.db')
    make_db('.OFprot', 'protDB_OF.db')

    os.chdir('..')