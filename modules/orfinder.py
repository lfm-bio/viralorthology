'''
This script launches orffinder with every genome (DNA output), sorts the seqs by length (from longer to shorter), gives an ID to every seq, 
translates them to protein and cleans the output (removes proteins that are annotated)
'''
import os
from Bio import SeqIO
from tqdm import tqdm
from modules import commands
from modules.misc import get_file_list
from modules.misc import sort_bioseqs_b_to_s
from modules.misc import get_bioseqs

def launch_orfinder(params):
    genomes = get_file_list()
    params = '-outfmt 1 ' + params # outfmt = DNA
    print('Running ORFfinder')
    for genome in tqdm(genomes):
        file_in, file_out = genome, genome.replace('.fasta', '.OForf')
        commands.run_orfinder(file_in, file_out, params)

def sort_genes_by_len():
    '''
    sorts seqs from bigger to smaller
    '''
    orfeomes = get_file_list('.OForf')
    for orfeome in orfeomes:
        genes = get_bioseqs(orfeome)
        genes = sort_bioseqs_b_to_s(genes)
        tmp_name = orfeome + '_tmp'
        SeqIO.write(genes, tmp_name, 'fasta-2line')
        os.remove(orfeome)
        os.rename(tmp_name, orfeome)

def get_ID(seq, n):
    ID = f'ORFINDER{n} '
    ID = seq.id[:seq.id.find(':')].replace('lcl|', ID)
    return ID

def edit_orfeomes():
    '''
    gives ID and format to every seq
    '''
    orfeomes = get_file_list('.OForf')
    n = 0
    for orfeome in orfeomes:
        tmp_name = orfeome + '_tmp'
        output = open(tmp_name, 'w', encoding='utf-8')
        for seq in SeqIO.parse(orfeome, 'fasta'):
            n += 1
            ID = get_ID(seq, n)
            position = seq.id[seq.id.find(':')+1:]
            position = ('..').join(position.split('-'))
            if position.startswith('c'):
                position = f'complement({position[1:]})'
            output.write(f'>{ID} [location={position}]\n{seq.seq}\n')
        output.close()
        os.remove(orfeome)
        os.rename(tmp_name, orfeome)

def translate():
    orfeomes = get_file_list('.OForf')
    for orfeome in orfeomes:
        prot_fasta = orfeome.replace('.OForf', '.OFprot')
        with open(prot_fasta, 'w', encoding='utf-8') as output:
            for seq in SeqIO.parse(orfeome, 'fasta'):
                output.write(f'>{seq.description}\n{seq.seq.translate(stop_symbol="")}\n')

def move_files():
    files = get_file_list('.OForf')
    for xfile in files:
        os.replace(xfile, f'../orfeomes/{xfile}')
    files = get_file_list('.OFprot')
    for xfile in files:
        os.replace(xfile, f'../proteomes/{xfile}')

def delete_seqs(OFprots_to_remove):
    of_files = get_file_list('.OFprot')
    for xfile in of_files:
        tmp_name = xfile + '_tmp'
        output = open(tmp_name, 'w', encoding='utf-8')
        for seq in SeqIO.parse(xfile, 'fasta'):
            if seq.id not in OFprots_to_remove:
                output.write(seq.format('fasta-2line'))
        output.close()
        os.remove(xfile)
        os.rename(tmp_name, xfile)

def clean_OFfiles():
    '''
    removes from OF output prots that are already annotated (in, ==, out)
    '''
    discarted = open('discarted.OF', 'w', encoding='utf-8')
    print('Cleaning ORFfinder output...')
    OFprots_to_remove = []
    of_files = get_file_list('.OFprot')
    for of_file in tqdm(of_files):
        anotated_fasta = of_file.replace('.OFprot', '.fasta')
        for seqOF in SeqIO.parse(of_file, 'fasta'):
            for seqAN in SeqIO.parse(anotated_fasta, 'fasta'):
                if seqOF.seq in seqAN.seq or seqOF.seq == seqAN.seq or seqAN.seq in seqOF.seq:
                    OFprots_to_remove.append(seqOF.id)
                    discarted.write(seqOF.format('fasta-2line'))
                    break
    discarted.close()
    delete_seqs(OFprots_to_remove)

def make_protDB_ingroup():
    fastas = get_file_list()
    protDB_ingroup = open('protDB_ingroup.fasta', 'w', encoding='utf-8')
    for group in fastas:
        for seq in SeqIO.parse(group, 'fasta'):
            protDB_ingroup.write(seq.format('fasta-2line'))
    protDB_ingroup.close()
    commands.makeblastdb_prot('protDB_ingroup.fasta')

def make_OF_db():
    with open('../orfeomes_OF.fasta', 'a', encoding='utf-8') as OF_db:
        OF_orfs = get_file_list('.OForf')
        for OFfile in OF_orfs:
            for seq in SeqIO.parse(OFfile, 'fasta'):
                OF_db.write(seq.format('fasta-2line'))

def main(params = '-ml 90 -s 0'):
    os.chdir('genomes')

    launch_orfinder(params)
    sort_genes_by_len()
    edit_orfeomes()
    make_OF_db()
    translate()
    move_files()
    os.chdir('../proteomes')
    clean_OFfiles()

    os.chdir('..')