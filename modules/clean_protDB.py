'''
This script removes from protDB all the seqs that are in some ortology group
'''
import os
from Bio import SeqIO
from modules.misc import get_file_list

def get_prots_in_group():
    '''
    OUT: list with the IDs of grouped proteins
    '''
    ortology_groups = get_file_list()
    prots_in_group = []
    for group in ortology_groups:
        for seq in SeqIO.parse(group, 'fasta'):
            prots_in_group.append(seq.id)
    return prots_in_group

def cleanDB():
    prots_in_group = get_prots_in_group()
    db_name = '../protDB.db'
    tmp_name = db_name + '_tmp'
    new_db = open(tmp_name, 'w')
    for seq in SeqIO.parse(db_name, 'fasta'):
        if seq.id not in prots_in_group:
            new_db.write(seq.format('fasta-2line'))
    new_db.close()
    os.remove(db_name)
    os.rename(tmp_name, db_name)

def main():

    os.chdir('orthology_groups')

    cleanDB()

    os.chdir('..')