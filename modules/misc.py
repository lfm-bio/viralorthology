import os
import multiprocessing
import time
from Bio import SeqIO

def combine_fastas(fasta_list, new_fasta_name, format = 'fasta-2line'):
    bioseqs = []
    for fasta in fasta_list:
        bioseqs += get_bioseqs(fasta)
    delete_files(fasta_list)
    SeqIO.write(bioseqs, new_fasta_name, format)

def delete_files(file_list):
    for xfile in file_list:
        os.remove(xfile)

def write_error_log(stage, file):
    error_log = open('../errors.log', 'a', encoding='utf-8')
    error_log.write(f'{stage} - {file}\n')
    error_log.close()

def get_first_bioseq_fasta(fasta):
    first_seq = list(SeqIO.parse(fasta, 'fasta'))[0]
    return first_seq

def delete_tmp_files(extensions):
    '''
    IN: list with extensions
    '''
    for file_ext in extensions:
        tmp_files = get_file_list(file_ext)
        for tmp_file in tmp_files:
            os.remove(tmp_file)

def check_ortology_group(query_genome, fasta): # this is ugly
    '''
    checks if the group already has a protein from that genome
    OUT: True if the group doesnt have a protein from that genome
    '''
    fasta_genomes = get_genomes_fasta(fasta)
    if query_genome in fasta_genomes:
        return False
    else:
        return True

def add_bioseq_to_fasta(bioseq, fasta):
    fasta = open(fasta, 'a', encoding='utf-8')
    fasta.write(bioseq.format('fasta-2line'))
    fasta.close()

def run_blastp(query_fasta, db_fasta, params = ''):
    n_cpus = get_cpu_number()
    err = os.system(f'blastp -query {query_fasta} -db {db_fasta} -out blastp.results -num_threads {n_cpus} {params}') #SALIDA BLASTP CAMBIO 
    if err:
        write_error_log('blastp', query_fasta)
        print(f'Error with {query_fasta}')
        time.sleep(3)

def get_blastp_hits():
    '''
    reads blast report
    OUT: list with the groups that gave a hit (from best to worst hit)
    [(evalue, hit)]
    '''
    read = False
    hits = []
    results = open('blastp.results', encoding='utf-8')
    for line in results:
        if line.find('No hits found') != -1: #no hit
            break
        elif line.startswith('Sequences producing significant alignments:'):
            read = True
            continue
        if read:
            if not line.strip():
                continue
            elif line.startswith('>'):
                break
            line = line.strip().split()
            fasta, evalue = line[0], float(line[-1]) #first: protID or fasta, last: evalue
            hit = (evalue, fasta)
            hits.append(hit)
    results.close()
    os.remove('blastp.results')
    return hits

def get_genome_bioseq(seq):
    genome = seq.description.split()[1]
    return genome

def makeblastdb_prot(fasta):
    err = os.system(f'makeblastdb -dbtype prot -in {fasta}')
    if err:
        print(f'Error with {fasta}')
        write_error_log('makeblastdb', fasta)
        time.sleep(3)

def align_fasta_muscle(fasta):
    '''
    OUT: aligned fasta file name
    '''
    aligned_fasta = fasta.replace('.fasta', '.muscle')
    if not os.path.isfile(aligned_fasta):
        os.system(f'muscle -align {fasta} -output {aligned_fasta}')
    return aligned_fasta

def make_hmm_hammer(aligned_fasta):
    '''
    OUT: hmm file name
    '''
    HMM = aligned_fasta.replace('.muscle', '.hmmbuild')
    os.system(f'hmmbuild {HMM} {aligned_fasta}')
    return HMM

def search_with_hmm(HMM, db_path, params):
    '''
    OUT: hmm search report file name
    '''
    search_hmm = HMM.replace('.hmmbuild', '.hmmsearch')
    os.system(f'hmmsearch -o {search_hmm} {params} {HMM} {db_path}')
    return search_hmm

def get_bioseqs(fasta):
    bioseqs = [seq for seq in SeqIO.parse(fasta, 'fasta')]
    return bioseqs

def sort_bioseqs_b_to_s(bioseqs):
    bioseqs.sort(key=lambda seq: len(seq.seq), reverse=True)
    return bioseqs

def run_orfinder(file_in, file_out, params):
    os.system(f'ORFfinder -in {file_in} -out {file_out} {params}')

def get_file_list(file_ext = '.fasta'):
    files = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith(file_ext)]
    files.sort()
    return files

def clean_protDB(prots_to_remove, protDB_path):
    '''
    removes genes added to a ortology group from protDB
    '''
    tmp_name = protDB_path + '_tmp'
    new_protDB = open(tmp_name, 'w', encoding='utf-8')

    for gene in SeqIO.parse(protDB_path, 'fasta'):
        if gene.id not in prots_to_remove:
            new_protDB.write(gene.format('fasta-2line'))

    new_protDB.close()
    os.remove(protDB_path)
    os.rename(tmp_name, protDB_path)

def get_nseqs(fasta):
    '''
    OUT: number of seqs in fasta file
    '''
    fasta = open(fasta, encoding='utf-8')
    n = 0
    for line in fasta:
        if line.startswith('>'):
            n += 1
    fasta.close()
    return n

def get_proteinid_genomeid(fasta):
    '''
    OUT: [(protID, genomeID), (protID, genomeID), ...]
    '''
    ids = []
    fasta = open(fasta, encoding='utf-8')
    for line in fasta:
        if line.startswith('>'):
            line = line[1:].split()
            protID, genomeID = line[0], line[1]
            ids.append((protID, genomeID))
    fasta.close()
    return ids

def get_protein_ids(fasta):
    '''
    OUT: list with the id of every protein in that fasta
    '''
    fasta = open(fasta, encoding='utf-8')
    ids = []
    for line in fasta:
        if line.startswith('>'):
            ids.append(line[1:].split()[0])
    fasta.close()
    return ids

def get_n_genomes(from_genomesfasta = False):
    '''
    OUT: number of genomes
    if from_genomesfasta it counts the number directly from genomes.fasta
    '''
    if from_genomesfasta or os.path.isfile('../first-round-protDB.db'):
        n_genomes =  get_nseqs('../genomes.fasta')
        return n_genomes
    n_genomes = len([xfile for xfile in os.listdir('../genomes') if xfile.endswith('.fasta')])
    if os.path.isfile('../genomes_old.fasta'): #Total number of genomes for update
        n_genomes += get_nseqs('../genomes_old.fasta')

    return n_genomes

def get_genomes_fasta(fasta):
    '''
    OUT: list with the IDs of the genomes which have a protein in the group
    '''
    fasta = open(fasta, encoding='utf-8')
    genomes_fasta = []
    for line in fasta:
        if line.startswith('>'):
            genome = line.split()[1]
            genomes_fasta.append(genome)
    fasta.close()
    return genomes_fasta

def get_cpu_number():
    return multiprocessing.cpu_count()

def get_ordered_files(file_list):
    '''
    OUT: list with fastas from biggest to smallest
    '''
    ordered_files = []
    for xfile in file_list:
        n_seqs = get_nseqs(xfile)
        ordered_files.append((n_seqs, xfile))
    ordered_files.sort(reverse = True)
    ordered_list = [xfile[1] for xfile in ordered_files]
    return ordered_list

def make_nseq_report(stage):
    output_path = 'n_seqs_groups.csv'

    if not os.path.isfile(output_path):
        write_header = True
    else:
        write_header = False
    output = open(output_path, 'a', encoding='utf-8')
    if write_header:
        output.write('file,stage,n_seqs\n')

    os.chdir('orthology_groups')
    fastas = get_file_list()
    for fasta in fastas:
        n_seqs = get_nseqs(fasta)
        output.write(f'{fasta},{stage},{n_seqs}\n')
    output.close()
    os.chdir('..')

def fasta_to_fasta2line(fasta):
    tmp = open('tmp', 'w', encoding='utf-8')
    for seq in SeqIO.parse(fasta, 'fasta'):
        tmp.write(seq.format('fasta-2line'))
    tmp.close()
    os.remove(fasta)
    os.rename('tmp', fasta)