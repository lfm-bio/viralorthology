#!/usr/bin/python3
import os
import sys
from Bio import SeqIO
import modules.split_data as split_data
import modules.orfinder as orfinder
import modules.make_protDB as make_protDB
import modules.search_paralogs as paralogs
import modules.blastp as blastp
import modules.HMM as HMM
import modules.synteny as synteny
from modules.misc import make_nseq_report
from modules.misc import combine_fastas

def split_params(params):
    for n, w  in enumerate(params):
        if w == '-query':
            query = params[n+1]
        elif w == '-mindate':
            mindate = params[n+1]
        elif w == '-minlenght':
            minlen = int(params[n+1])
    return query, mindate, minlen

def download_seqs(query, mindate):
    err = 1
    while err != 0:
        err = os.system(f'esearch -db nuccore -query "{query} AND viruses[filter]" | efilter -mindate {mindate} | efetch -format fasta > new_genomes.fasta')

def filter_by_len(minlen):
    tmp = open('tmp.fasta', 'w', encoding='utf-8')
    for seq in SeqIO.parse('new_genomes.fasta', 'fasta'):
        if 'partial genome' in seq.description:
            continue
        if len(seq.seq) >= minlen:
            tmp.write(seq.format('fasta'))
    tmp.close()
    os.remove('new_genomes.fasta')
    os.rename('tmp.fasta', 'new_genomes.fasta')

def remove_identicals():
    '''
    removes genomes with same accs number or same sequence
    '''
    ids = []
    seqs = []
    for seq in SeqIO.parse('genomes_old.fasta', 'fasta'):
        ids.append(seq.id)
        seqs.append(seq.seq)
    tmp = open('tmp.fasta', 'w', encoding='utf-8')
    for seq in SeqIO.parse('new_genomes.fasta', 'fasta'):
        if seq.id not in ids and seq.seq not in seqs:
            tmp.write(seq.format('fasta'))
    tmp.close()
    os.remove('new_genomes.fasta')
    os.rename('tmp.fasta', 'new_genomes.fasta')

def remove_wrong_genomes():
    '''
    prints description of new seqs for the user to choose which to remove
    '''
    ids = []
    print('N - Seq. Description')
    for n, seq in enumerate(SeqIO.parse('new_genomes.fasta', 'fasta')):
        ids.append(seq.id)
        description = (' ').join(seq.description.split()[1:])
        print(f'{n} - {description}')
    index_to_remove = input('Enter N of seqs to remove separated by spaces: ').split()
    index_to_remove = [int(N) for N in index_to_remove]
    ids_to_remove = [seqid for n, seqid in enumerate(ids) if n in index_to_remove]
    with open('tmp.fasta', 'w', encoding='utf-8') as tmp:
        for seq in SeqIO.parse('new_genomes.fasta', 'fasta'):
            if seq.id not in ids_to_remove:
                tmp.write(seq.format('fasta'))
    os.remove('new_genomes.fasta')
    os.rename('tmp.fasta', 'genomes.fasta')

def download_proteome_orfeome():
    ids = []
    orfeomes = open('orfeomes.fasta', 'w', encoding='utf-8')
    proteomes = open('proteomes.fasta', 'w', encoding='utf-8')
    for seq in SeqIO.parse('genomes.fasta', 'fasta'):
        ids.append(seq.id)
    for seqid in ids:
        while True:
            os.system(f'efetch -db nuccore -id "{seqid}" -format fasta_cds_aa > new_proteome.fasta') #proteome
            os.system(f'efetch -db nuccore -id "{seqid}" -format fasta_cds_na > new_orfeome.fasta') #orfeome
            if os.path.exists('new_proteome.fasta') and os.path.exists('new_orfeome.fasta'):
                break
            print(f'Trying again with {seqid}')
        for seq in SeqIO.parse('new_proteome.fasta', 'fasta'):
            proteomes.write(seq.format('fasta'))
        for seq in SeqIO.parse('new_orfeome.fasta', 'fasta'):
            orfeomes.write(seq.format('fasta'))
        os.remove('new_orfeome.fasta')
        os.remove('new_proteome.fasta')
    orfeomes.close()
    proteomes.close()

def pipeline():
    params_HMMsearch = False
    split_data.main()

    orfinder.main()
    paralogs.main()
    make_protDB.main()

    HMM.main(params_HMMsearch)
    make_nseq_report('2-HMM')

    blastp.main('protDB.db')
    blastp.main('protDB_OF.db')
    make_nseq_report('3-Blastp')

    HMM.main(params_HMMsearch, True)
    make_nseq_report('4-HMM')

    synteny.main()

def main(params = ''):
    if os.path.exists('synteny.txt'):
        os.remove('synteny.txt')
    os.rename('genomes.fasta', 'genomes_old.fasta')
    os.rename('proteomes.fasta', 'proteomes_old.fasta')
    os.rename('orfeomes.fasta', 'orfeomes_old.fasta')
    os.rename('protDB.db', 'protDB_old.db')

    #params = '-query baculovirus -mindate 2022 -minlenght 50000'
    params = sys.argv[2:]
    query, mindate, minlen = split_params(params)
    print('Downloading new genomes...')
    download_seqs(query, mindate)
    filter_by_len(minlen)
    remove_identicals()
    remove_wrong_genomes()
    print('Downloading new proteomes and orfeomes...')
    download_proteome_orfeome()
    pipeline()

    combine_fastas(['genomes_old.fasta', 'genomes.fasta'], 'genomes.fasta', 'fasta')
    combine_fastas(['proteomes_old.fasta', 'proteomes.fasta'], 'proteomes.fasta', 'fasta')
    combine_fastas(['orfeomes_old.fasta', 'orfeomes.fasta'], 'orfeomes.fasta', 'fasta')
    combine_fastas(['protDB_old.db', 'protDB.db'], 'protDB.db')
    os.system('rm -r genomes')
    os.system('rm -r orfeomes')
    os.system('rm -r proteomes')
    os.remove('n_seqs_groups.csv')
    os.remove('protDB_OF.db')