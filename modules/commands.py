import os
import subprocess
import sys
import time
from modules.misc import get_cpu_number
from modules.misc import write_error_log

DONT_PRINT_OUTPUT = '> /dev/null'

def make_blastdb_nt(fasta):
    err = os.system(f'makeblastdb -dbtype nucl -in {fasta} {DONT_PRINT_OUTPUT}')
    if err:
        sys.exit(1)

def makeblastdb_prot(fasta):
    err = os.system(f'makeblastdb -dbtype prot -in {fasta} {DONT_PRINT_OUTPUT}')
    if err:
        print(f'Error with {fasta}')
        write_error_log('makeblastdb', fasta)
        time.sleep(3)

def blast_filter_genomes(f1, f2):
    os.system(f"blastn -query {f1} -db {f2} -out blast.results -outfmt '7 delim=@ qcovs pident' {DONT_PRINT_OUTPUT}")

def blastp(query_fasta, db_fasta, params = ''):
    n_cpus = get_cpu_number()
    err = os.system(f'blastp -query {query_fasta} -db {db_fasta} -out blastp.results -num_threads {n_cpus} {params} {DONT_PRINT_OUTPUT}')
    if err:
        write_error_log('blastp', query_fasta)
        print(f'Error with {query_fasta}')
        time.sleep(3)

def download_genome(seqid):
    return os.system(f'efetch -db nuccore -id "{seqid}" -format fasta > {seqid}.genome')

def download_proteome(seqid):
    return os.system(f'efetch -db nuccore -id "{seqid}" -format fasta_cds_aa > {seqid}.proteome')

def download_orfeome(seqid):
    return os.system(f'efetch -db nuccore -id "{seqid}" -format fasta_cds_na > {seqid}.orfeome')

def proteinortho(proteomes, params):
    with open(os.devnull, 'wb') as devnull:
        subprocess.check_call(f'proteinortho {proteomes} {params}'.split(), stdout=devnull, stderr=subprocess.STDOUT)
    with open(os.devnull, 'wb') as devnull:
        subprocess.check_call(f'po_grab_proteins -tofiles myproject.proteinortho.tsv -exact {proteomes}'.split(), stdout=devnull, stderr=subprocess.STDOUT)
    # print(proteomes)
    # quit()
    # os.system(f'proteinortho {proteomes} {params}')
    # os.system(f'po_grab_proteins -tofiles myproject.proteinortho.tsv -exact {proteomes}')

def align_muscle(fasta, aligned_fasta):
    with open(os.devnull, 'wb') as devnull:
        subprocess.check_call(f'muscle -align {fasta} -output {aligned_fasta}'.split(), stdout=devnull, stderr=subprocess.STDOUT)

def build_hmm(HMM, aligned_fasta):
    os.system(f'hmmbuild {HMM} {aligned_fasta} {DONT_PRINT_OUTPUT}')

def search_hmm(search_hmm, params, HMM, db_path):
    os.system(f'hmmsearch -o {search_hmm} {params} {HMM} {db_path} {DONT_PRINT_OUTPUT}')

def run_orfinder(file_in, file_out, params):
    os.system(f'ORFfinder -in {file_in} -out {file_out} {params} {DONT_PRINT_OUTPUT}')

def fasta_to_a2m(output_muscle, output_a2m):
    with open(os.devnull, 'wb') as devnull:
        subprocess.check_call(f'reformat_msa fas a2m {output_muscle} {output_a2m}'.split(), stdout=devnull, stderr=subprocess.STDOUT)

def hmm_make(output_a2m, n_genomes):
    with open(os.devnull, 'wb') as devnull:
        subprocess.check_call(f'hhmake -i {output_a2m} -id 100 -diff {n_genomes}'.split(), stdout=devnull, stderr=subprocess.STDOUT)

def hmm_align(g1, g2, report):
    with open(os.devnull, 'wb') as devnull:
        subprocess.check_call(f'hhalign -i {g1} -t {g2} -o {report} -glob'.split(), stdout=devnull, stderr=subprocess.STDOUT)