'''
renames the orthology groups going from bigger to smaller
'''
import os
from modules.misc import get_file_list
from modules.misc import get_ordered_files

def tmp_names():
    '''
    gives a temporal name to all fasta files (to avoid overwriting files)
    '''
    files = get_file_list()
    tmp_files = []
    for old_name in files:
        tmp_name = old_name + '_tmp'
        tmp_files.append(tmp_name)
        os.rename(old_name, tmp_name)
    return tmp_files

def rename_groups():
    file_list = tmp_names()
    ordered_files = get_ordered_files(file_list)

    n_figures = len(str(len(file_list)))
    for n, fasta in enumerate(ordered_files):
        old_name = fasta
        n = str(n)
        while len(n) < n_figures:
            n = '0' + n
        new_name = f'{n}.fasta'
        os.rename(old_name, new_name)

def main():
    os.chdir('orthology_groups')

    rename_groups()

    os.chdir('..')