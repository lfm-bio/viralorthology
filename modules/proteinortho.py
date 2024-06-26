import os
from modules import commands
from modules.misc import get_file_list
from modules.misc import delete_tmp_files

def launch_proteinortho(params):
    proteomes = get_file_list()
    proteomes = (' ').join(proteomes)
    print('Running proteinortho...')
    commands.proteinortho(proteomes, params)

def move_files():
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

def main(params):

    if not params:
        params = '--p=blastp+'

    os.chdir('proteomes')

    launch_proteinortho(params)
    move_files()
    clean_tmpfiles()

    os.chdir('..')
