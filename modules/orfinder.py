import os
from Bio import SeqIO
from tqdm import tqdm
from modules import commands
from modules.misc import get_file_list

def launch_orfinder(params):
    genomes = get_file_list()
    params = '-outfmt 1 ' + params # outfmt = DNA
    print('Running ORFfinder...')
    for genome in tqdm(genomes):
        file_in, file_out = genome, genome.replace('.fasta', '.OForf')
        commands.run_orfinder(file_in, file_out, params)

def sort_genes_by_len_big_to_small():
    orfeomes = get_file_list('.OForf')
    for orfeome in orfeomes:
        genes = list(SeqIO.parse(orfeome, 'fasta'))
        genes.sort(key=lambda seq: len(seq.seq), reverse=True)
        SeqIO.write(genes, 'tmp', 'fasta')
        os.remove(orfeome)
        os.rename('tmp', orfeome)

def give_id_format_seqs():

    def get_id(seq, n):
        id_ = f'ORFINDER{n} '
        id_ = seq.id[:seq.id.find(':')].replace('lcl|', id_)
        return id_

    orfeomes = get_file_list('.OForf')
    n = 0
    for orfeome in orfeomes:
        with open('tmp', 'w', encoding='utf-8') as output:
            for seq in SeqIO.parse(orfeome, 'fasta'):
                n += 1
                id_ = get_id(seq, n)
                position = seq.id[seq.id.find(':')+1:]
                position = ('..').join(position.split('-'))
                if position.startswith('c'):
                    position = f'complement({position[1:]})'
                output.write(f'>{id_} [location={position}]\n{seq.seq}\n')
        os.remove(orfeome)
        os.rename('tmp', orfeome)

def translate_orfeomes():
    fastas = get_file_list('.OForf')
    for fasta in fastas:
        prot_fasta = fasta.replace('.OForf', '.OFprot')
        with open(prot_fasta, 'w', encoding='utf-8') as output:
            for seq in SeqIO.parse(fasta, 'fasta'):
                output.write(f'>{seq.description}\n{seq.seq.translate(stop_symbol="")}\n')

def move_fastas():
    fastas = get_file_list('.OForf')
    for fasta in fastas:
        os.replace(fasta, f'../orfeomes/{fasta}')
    fastas = get_file_list('.OFprot')
    for fasta in fastas:
        os.replace(fasta, f'../proteomes/{fasta}')

def clean_of_proteomes():
    '''
    removes from OF output prots that are already annotated (in, ==, out)
    '''

    def remove_of_prots(prots_to_remove):
        fastas = get_file_list('.OFprot')
        for fasta in fastas:
            with open('tmp', 'w', encoding='utf-8') as output:
                for seq in SeqIO.parse(fasta, 'fasta'):
                    if seq.id not in prots_to_remove:
                        output.write(seq.format('fasta-2line'))
            os.remove(fasta)
            os.rename('tmp', fasta)

    of_prots_to_remove = []
    of_files = get_file_list('.OFprot')
    discarted = open('discarted.OF', 'w', encoding='utf-8')
    print('Cleaning ORFfinder output...')
    for of_file in tqdm(of_files):
        annotated_fasta = of_file.replace('.OFprot', '.fasta')
        for of_seq in SeqIO.parse(of_file, 'fasta'):
            for ann_seq in SeqIO.parse(annotated_fasta, 'fasta'):
                if of_seq.seq in ann_seq.seq or of_seq.seq == ann_seq.seq or ann_seq.seq in of_seq.seq:
                    of_prots_to_remove.append(of_seq.id)
                    discarted.write(of_seq.format('fasta'))
                    break
    discarted.close()
    remove_of_prots(of_prots_to_remove)

def make_of_db():
    with open('../orfeomes_OF.fasta', 'a', encoding='utf-8') as of_db:
        fastas = get_file_list('.OForf')
        for fasta in fastas:
            for seq in SeqIO.parse(fasta, 'fasta'):
                of_db.write(seq.format('fasta'))

def main(params):
    if not params:
        params = '-ml 90 -s 0'

    os.chdir('genomes')

    launch_orfinder(params)
    sort_genes_by_len_big_to_small()
    give_id_format_seqs()
    make_of_db()
    translate_orfeomes()
    move_fastas()
    os.chdir('../proteomes')
    clean_of_proteomes()

    os.chdir('..')
