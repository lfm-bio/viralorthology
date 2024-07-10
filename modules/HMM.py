import os
from Bio import SeqIO
from tqdm import tqdm
from modules.misc import remove_from_prot_db
from modules.misc import get_genomes_fasta
from modules.misc import get_n_genomes
from modules.misc import get_nseqs
from modules.misc import get_ordered_files
from modules.misc import get_file_list
from modules.misc import align_fasta_muscle
from modules.misc import make_hmm_hammer
from modules.misc import search_with_hmm
from modules.misc import delete_tmp_files

def get_hits_hmmsearch(fasta):
    '''
    OUT: dict with the genes gave a hit (only the best hit per genome)
    dict[genome] = gene
    '''
    genes = {}
    hmmsearch_report = open(fasta.replace('.fasta', '.hmmsearch'), encoding='utf-8')
    read = False
    for line in hmmsearch_report:
        if line.strip().startswith('E-value  score  bias'):
            read = True
            continue
        if read:
            if line.strip().startswith('------ inclusion threshold ------') or line.strip().startswith('Domain annotation for each sequence') or not line.strip():
                break
            if line.strip().startswith('-'):
                continue
            line = line.strip().split()
            gene, genome = line[8], line[9]
            if genome not in genes: #i get only the best hit per genome
                genes[genome] = gene
    hmmsearch_report.close()
    return genes

def get_new_genes(fasta):
    '''
    OUT: list of genes to add to the orthology group
    checks theres not another gene from the same genome in that orthology group
    '''


    hits = get_hits_hmmsearch(fasta) #dict[genome] = gene
    genomes_fasta = get_genomes_fasta(fasta)
    new_genes = [hit for genome, hit in hits.items() if genome not in genomes_fasta]

    return new_genes

def add_new_genes_to_fasta(fasta, new_genes):
    with open(fasta, 'a', encoding='utf-8') as fasta_handle:
        for gene in SeqIO.parse('../protDB.db', 'fasta'):
            if gene.id in new_genes:
                fasta_handle.write(gene.format('fasta'))

def align_build_search(fasta, search_params = ''):
    align_fasta_muscle(fasta)
    make_hmm_hammer(fasta)
    search_with_hmm(fasta, '../protDB.db', search_params)

def hmm(search_params):
    prot_db_path = '../protDB.db'
    fastas = get_file_list()
    fastas = get_ordered_files(fastas)
    n_genomes = get_n_genomes()
    print('Searching for new genes with HMM...')
    for fasta in tqdm(fastas):
        while True:
            n_genes = get_nseqs(fasta)
            if n_genes == n_genomes: #if the group already has one protein per genome
                break
            align_build_search(fasta, search_params)
            new_genes = get_new_genes(fasta)
            if new_genes:
                add_new_genes_to_fasta(fasta, new_genes)
                remove_from_prot_db(new_genes, prot_db_path)
                os.remove(fasta.replace('.fasta', '.muscle')) #so it has to align the fasta again the next iteration
            else: #no new proteins found
                break

def main(search_params):
    os.chdir('orthology_groups')

    hmm(search_params)
    delete_tmp_files(['.hmmbuild', '.hmmsearch', '.muscle'])

    os.chdir('..')
