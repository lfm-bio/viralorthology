import os
from Bio import SeqIO
from tqdm import tqdm
from modules.HMM import get_hits_hmmsearch
from modules.HMM import align_build_search
from modules.misc import delete_tmp_files
from modules.misc import get_file_list
from modules.misc import get_ordered_files
from modules.misc import get_nseqs

def remake_group(fasta, hits):
    with open('tmp', 'w', encoding='utf-8') as tmp_file:
        for seq in SeqIO.parse(fasta, 'fasta'):
            if seq.id in hits:
                tmp_file.write(seq.format('fasta'))

    os.remove(fasta)
    if get_nseqs('tmp') < 2: #group with only one seq
        os.remove('tmp')
    else:
        os.rename('tmp', fasta)

def clean_groups():
    fastas = get_file_list()
    fastas = get_ordered_files(fastas)
    print('Cleaning groups...')
    for fasta in tqdm(fastas):
        intial_nseqs = get_nseqs(fasta)
        align_build_search(fasta)
        hits = list(get_hits_hmmsearch(fasta).values())
        remake_group(fasta, hits)
        final_nseqs = get_nseqs(fasta)
        if intial_nseqs != final_nseqs: #the groups has been modified, align again to hmmvshmm
            os.remove(fasta.replace('.fasta', '.muscle'))

def main():
    os.chdir('orthology_groups')

    clean_groups()
    delete_tmp_files(['.hmmbuild', '.hmmsearch'])

    os.chdir('..')