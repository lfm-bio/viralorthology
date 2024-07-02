import os
from Bio import SeqIO
from modules.misc import get_file_list

def get_prots_in_group():
    '''
    OUT: list with the IDs of grouped proteins
    '''
    fastas = get_file_list()
    prots_in_group = []
    for fasta in fastas:
        for seq in SeqIO.parse(fasta, 'fasta'):
            prots_in_group.append(seq.id)
    return prots_in_group

def cleanDB():
    prots_in_group = get_prots_in_group()
    db_name = '../protDB.db'
    with open('../tmp', 'w', encoding='utf-8') as new_db:
        for seq in SeqIO.parse(db_name, 'fasta'):
            if seq.id not in prots_in_group:
                new_db.write(seq.format('fasta'))
    os.remove(db_name)
    os.rename('../tmp', db_name)

def main():

    os.chdir('orthology_groups')

    cleanDB()

    os.chdir('..')
