'''
This script launches proteinortho and grab_proteins
'''
import os
from modules.misc import get_file_list
from modules.misc import delete_tmp_files

def proteinortho(params):
    '''
    launches proteinortho and grab_proteins
    '''
    proteomes = get_file_list()
    proteomes = (' ').join(proteomes)
    os.system(f'proteinortho {proteomes} {params}')
    os.system(f'po_grab_proteins -tofiles myproject.proteinortho.tsv -exact {proteomes}')

def move_files():
    '''
    moves orthology groups to orthology_groups folder
    '''
    os.mkdir('../orthology_groups')
    groups = [xfile for xfile in os.listdir(os.curdir) if xfile.find('OrthoGroup') != -1]
    for group in groups:
        os.replace(group, f'../orthology_groups/{group}')

def clean_tmpfiles():
    tmpfiles = [xfile for xfile in os.listdir(os.curdir) if xfile.find('myproject') != -1 or xfile.endswith('.dmnd')]
    for tmpfile in tmpfiles:
        os.remove(tmpfile)
    os.chdir('../proteomes')
    delete_tmp_files(['.phr', '.pin', '.psq'])

def main(params = '--p=blastp+'):

    os.chdir('proteomes')

    proteinortho(params)
    move_files()
    clean_tmpfiles()

    os.chdir('..')