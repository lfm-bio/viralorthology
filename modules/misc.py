import os
import sys
import time
from datetime import timedelta, datetime
import multiprocessing
from Bio import SeqIO
from modules import commands
from classes.Args import Args

def get_args():

    def check_for_prohibited_params(args):
        '''
        checks if the user changed some prohibited parameter
        '''
        ok = True
        if args.proteinortho:
            if '--project' in args.proteinortho:
                ok = False
        if args.blastp:
            if '-query' in args.blastp or '-db' in args.blastp or '-out' in args.blastp:
                ok = False
        if args.orffinder:
            if '-in' in args.orffinder or '-out' in args.orffinder or '-outfmt' in args.orffinder:
                ok = False
        if not ok:
            print('Some prohibited parameter was changed')
            print('visit github for more information')
            sys.exit(1)

    args_dict = {}
    args = Args()

    software_list = [
        '--proteinortho',
        '--orffinder',
        '--blastp',
        '--HMMsearch'
        ]

    software = False
    for param in sys.argv[1:]:
        if param in software_list:
            software = param.strip('-')
            args_dict[software] = []
            continue
        if software:
            args_dict[software].append(param)

    for software, params in args_dict.items():
        setattr(args, software.strip('-'), (' ').join(params))

    check_for_prohibited_params(args)

    return args

def combine_fastas(fasta_list, new_fasta_name, file_format = 'fasta'):
    bioseqs = []
    for fasta in fasta_list:
        bioseqs += SeqIO.parse(fasta, 'fasta')
    delete_files(fasta_list)
    SeqIO.write(bioseqs, new_fasta_name, format=file_format)

def delete_files(file_list):
    for xfile in file_list:
        os.remove(xfile)

def write_error_log(stage, file):
    with open('../errors.log', 'a', encoding='utf-8') as error_log:
        error_log.write(f'{stage} - {file}\n')

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

def delete_tmp_files_final():
    os.chdir('orthology_groups')
    files = [xfile for xfile in os.listdir(os.curdir) if xfile.startswith('protDB')]
    for xfile in files:
        os.remove(xfile)
    os.chdir('..')

def check_ortology_group(query_genome, fasta): # this is ugly
    '''
    checks if the group already has a protein from that genome
    OUT: True if the group doesnt have a protein from that genome
    '''
    fasta_genomes = get_genomes_fasta(fasta)
    if query_genome in fasta_genomes:
        return False
    return True

def add_bioseq_to_fasta(bioseq, fasta):
    with open(fasta, 'a', encoding='utf-8') as fasta_handle:
        fasta_handle.write(bioseq.format('fasta'))

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
        if line.startswith('Sequences producing significant alignments:'):
            read = True
            continue
        if read:
            if not line.strip():
                continue
            if line.startswith('>'):
                break
            line = line.strip().split()
            fasta, evalue = line[0], float(line[-1]) #first: protID or fasta, last: evalue
            hit = (evalue, fasta)
            hits.append(hit)
    results.close()
    os.remove('blastp.results')
    return hits

def get_genome_id_from_bioseq(seq):
    genome = seq.description.split()[1]
    return genome

def align_fasta_muscle(fasta):
    aligned_fasta = fasta.replace('.fasta', '.muscle')
    if not os.path.isfile(aligned_fasta):
        commands.align_muscle(fasta, aligned_fasta)

def make_hmm_hammer(fasta):
    aligned_fasta = fasta.replace('.fasta', '.muscle')
    hmm = fasta.replace('.fasta', '.hmmbuild')
    commands.build_hmm(hmm, aligned_fasta)

def search_with_hmm(fasta, db_path, params):
    hmm = fasta.replace('.fasta', '.hmmbuild')
    search_hmm = fasta.replace('.fasta', '.hmmsearch')
    commands.search_hmm(search_hmm, params, hmm, db_path)

def get_file_list(file_ext = '.fasta'):
    files = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith(file_ext)]
    files.sort()
    return files

def get_all_genome_ids():
    genome_ids = []
    for seq in SeqIO.parse('../genomes.fasta', 'fasta'):
        genome_ids.append(seq.id[:seq.id.rfind('.')])
    return genome_ids

def remove_from_prot_db(prots_to_remove, prot_db_path):
    '''
    removes genes added to a orthology group from protDB
    '''
    with open('tmp', 'w', encoding='utf-8') as new_prot_db:
        for gene in SeqIO.parse(prot_db_path, 'fasta'):
            if gene.id not in prots_to_remove:
                new_prot_db.write(gene.format('fasta'))

    os.remove(prot_db_path)
    os.rename('tmp', prot_db_path)

def get_nseqs(fasta):
    '''
    OUT: number of seqs in fasta file
    '''
    if os.path.isfile(fasta):
        n_seqs = len(list(SeqIO.parse(fasta, 'fasta')))
        return n_seqs
    return 0

def get_proteinid_genomeid(fasta):
    '''
    OUT: [(protID, genomeID), (protID, genomeID), ...]
    '''
    ids = []
    with open(fasta, encoding='utf-8') as fasta:
        for line in fasta:
            if line.startswith('>'):
                line = line[1:].split()
                protID, genomeID = line[0], line[1]
                ids.append((protID, genomeID))
    return ids

def get_protein_ids(fasta):
    '''
    OUT: list with the id of every protein in that fasta
    '''
    ids = []
    with open(fasta, encoding='utf-8') as fasta:
        for line in fasta:
            if line.startswith('>'):
                ids.append(line[1:].split()[0])
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
    if os.path.isfile('../genomes_old.fasta'): # total number of genomes for update
        n_genomes += get_nseqs('../genomes_old.fasta')

    return n_genomes

def get_genomes_fasta(fasta):
    '''
    OUT: list with the IDs of the genomes which have a protein in the group
    '''
    genomes_fasta = []
    with open(fasta, encoding='utf-8') as fasta_handle:
        for line in fasta_handle:
            if line.startswith('>'):
                genome_id = line.split()[1]
                genomes_fasta.append(genome_id)
    return genomes_fasta

def get_cpu_number():
    return multiprocessing.cpu_count()

def get_ordered_files(fastas):
    '''
    OUT: list with fastas from biggest to smallest
    '''
    ordered_files = sorted([(get_nseqs(fasta), fasta) for fasta in fastas], reverse=True)
    ordered_list = [xfile[1] for xfile in ordered_files]
    return ordered_list

def make_nseq_report(stage):
    output_path = 'n_seqs_groups.csv'

    write_header = not bool(os.path.isfile(output_path))
    with open(output_path, 'a', encoding='utf-8') as output:
        if write_header:
            output.write('file,stage,n_seqs\n')

        os.chdir('orthology_groups')
        fastas = get_file_list()
        for fasta in fastas:
            n_seqs = get_nseqs(fasta)
            output.write(f'{fasta},{stage},{n_seqs}\n')
    os.chdir('..')

def fasta_to_fasta2line(fasta):
    with open('tmp', 'w', encoding='utf-8') as tmp:
        for seq in SeqIO.parse(fasta, 'fasta'):
            tmp.write(seq.format('fasta-2line'))
    tmp.close()
    os.remove(fasta)
    os.rename('tmp', fasta)

def write_log(args, start):
    end = time.time()
    elapsed_time = end - start
    elapsed_time = str(timedelta(seconds=elapsed_time))
    date = datetime.today().strftime('%Y-%m-%d')
    with open('log.txt', 'w', encoding="utf-8") as log:
        log.write(f'{date}\nElapsed time: {elapsed_time}\n')
        log.write('Parameters:\n')
        log.write(args.get_all_params())

def delete_final_files():
    os.system('rm -r genomes')
    os.system('rm -r orfeomes')
    os.system('rm -r proteomes')
    if os.path.exists('first-round-n_seqs_groups.csv'):
        os.remove('first-round-n_seqs_groups.csv')
    os.remove('protDB_OF.db')

def merge_protDBs(db):
    first_round_db_name = f'first-round-{db}'
    if not os.path.isfile(first_round_db_name): # no filtered genomes
        return False

    with open(db, 'a', encoding='utf-8') as new_db:
        with open(first_round_db_name, encoding='utf-8') as first_round_db:
            for line in first_round_db:
                new_db.write(line)
    os.remove(first_round_db_name)

def get_unique_gene_bioseq_by_id(gene_id):
    for seq in SeqIO.parse('../protDB.db', 'fasta'):
        if seq.id == gene_id:
            return seq

def check_files():
    if not os.path.isfile('genomes.fasta'):
        print('genomes.fasta not found, run viralorthology -download_seqs')
        sys.exit(1)
    if not os.path.isfile('orfeomes.fasta'):
        print('orfeomes.fasta not found, run viralorthology -download_seqs')
        sys.exit(1)
    if not os.path.isfile('proteomes.fasta'):
        print('proteomes.fasta not found, run viralorthology -download_seqs')
        sys.exit(1)

def check_filtered_genomes():

    def check_filtered():
        if os.path.isdir('filtered'):
            return True
        return False

    os.chdir('genomes')
    filtered = check_filtered()
    os.chdir('..')

    return filtered

def prepare_second_round():

    def delete_move_files():

        def move_files():
            os.chdir('filtered')
            files = list(os.listdir(os.curdir))
            for xfile in files:
                os.replace(xfile, f'../{xfile}')
            os.chdir('..')
            os.rmdir('filtered')

        for folder in ['genomes', 'orfeomes', 'proteomes']:
            os.chdir(folder)
            files = [xfile for xfile in os.listdir(os.curdir) if os.path.isfile(xfile)]
            for xfile in files:
                os.remove(xfile)
            move_files()
            os.chdir('..')

    delete_move_files()
    os.rename('protDB_OF.db', 'first-round-protDB_OF.db')
    os.rename('protDB.db', 'first-round-protDB.db')