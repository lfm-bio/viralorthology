import os
from Bio import SeqIO
from modules import commands
from modules.misc import remove_from_prot_db
from modules.misc import get_nseqs
from modules.misc import get_n_genomes
from modules.misc import get_ordered_files
from modules.misc import get_genome_id_from_bioseq
from modules.misc import get_blastp_hits
from modules.misc import get_file_list
from modules.misc import add_bioseq_to_fasta
from modules.misc import check_ortology_group
from modules.misc import delete_tmp_files

def make_db_protsingroup():
    '''
    makes protDB_ingroup.fasta and makeblastdb
    every seq has the description: >groupname.fasta for easier later use
    '''
    fastas = get_file_list()
    fastas = get_ordered_files(fastas)
    n_genomes = get_n_genomes()
    with open('protDB_ingroup.fasta', 'w', encoding='utf-8') as prot_db_ingroup:
        for fasta in fastas:
            n_genes = get_nseqs(fasta)
            if n_genes < n_genomes: #only adds prots from groups that are not full
                for seq in SeqIO.parse(fasta, 'fasta'):
                    prot_db_ingroup.write(f'>{fasta}\n{seq.seq}\n')
    commands.makeblastdb_prot('protDB_ingroup.fasta')

def run_read_blastp(params):
    '''
    launches blastp
    OUT: best hit / False
    '''
    commands.blastp('seq.query', 'protDB_ingroup.fasta', params)

    os.remove('seq.query')
    hits = get_blastp_hits()
    if hits:
        return hits[0] #best hit
    return False

def add_seqs(hits_all):
    prots_to_remove = []
    for genome in hits_all:
        for fasta in hits_all[genome]:
            _, seq = hits_all[genome][fasta]
            add_to_fasta = check_ortology_group(genome, fasta)
            if add_to_fasta:
                add_bioseq_to_fasta(seq, fasta)
                prots_to_remove.append(seq.id)
                if os.path.isfile(fasta.replace('.fasta', '.muscle')):
                    os.remove(fasta.replace('.fasta', '.muscle'))
    return prots_to_remove

def submain(params, prot_db_name):
    prot_db = f'../{prot_db_name}'
    prot_db_nseqs = get_nseqs(prot_db) #just to calculate the % of the process
    hits_all = {}
    for n, query in enumerate(SeqIO.parse(prot_db, 'fasta')):
        print(f'\r{round((n+1)/prot_db_nseqs*100, 1)}% - Blasting {query.id}', end='')
        query_genome = get_genome_id_from_bioseq(query)
        if query_genome not in hits_all:
            hits_all[query_genome] = {} #dict[genome] = dict[fasta] = (evalue, query)
        SeqIO.write(query, 'seq.query', 'fasta')
        best_hit = run_read_blastp(params)
        if best_hit:
            evalue, fasta = best_hit
            if fasta in hits_all[query_genome] and evalue > hits_all[query_genome][fasta][0]: #the last hit from that genome with that group was better (higher evalue)
                continue
            hits_all[query_genome][fasta] = (evalue, query)
    prots_to_remove = add_seqs(hits_all)
    remove_from_prot_db(prots_to_remove, prot_db)
    print() # final \n

def main(prot_db, params = ''):
    if '-word_size' not in params:
        params += ' -word_size 2'

    os.chdir('orthology_groups')

    make_db_protsingroup() #db with all the grouped genes
    submain(params, prot_db)
    delete_tmp_files(['.phr', '.pin', '.psq'])
    os.remove('protDB_ingroup.fasta')

    os.chdir('..')
