'''
This script removes the proteins that are not found using the HMM made with the orthology group
if the final group has less than 2 seqs, the group is deleted
'''
import os
from Bio import SeqIO
from modules.HMM import get_hits_hmmsearch
from modules.HMM import align_build_search
from modules.misc import delete_tmp_files
from modules.misc import get_file_list
from modules.misc import get_ordered_files

def remake_group(group, hits):
    '''
    OUT: True if the group has been modified (-> delete .muscle)
    else False: we will use the same .muscle in HMMvsHMM
    '''
    tmp_file = open('temp', 'w', encoding='utf-8')
    n_seqs = 0
    delete_muscle = False
    for seq in SeqIO.parse(group, 'fasta'):
        if seq.id in hits:
            tmp_file.write(seq.format('fasta-2line'))
            n_seqs += 1
        else: #delete seq from group -> align again to HMMvsHMM
            delete_muscle = True
    tmp_file.close()

    os.remove(group)
    if n_seqs < 2: #group with only one seq
        os.remove('temp')
    else:
        os.rename('temp', group)

    return delete_muscle

def clean_groups():
    groups = get_file_list()
    groups = get_ordered_files(groups)
    for group in groups:
        results = align_build_search(group)
        hits = list(get_hits_hmmsearch(results).values())
        groups_been_modified = remake_group(group, hits)
        if groups_been_modified: #the groups has been modified, align again to hmmvshmm
            os.remove(group.replace('.fasta', '.muscle'))

def main():
    os.chdir('orthology_groups')

    clean_groups()
    delete_tmp_files(['.hmmbuild', '.hmmsearch'])

    os.chdir('..')