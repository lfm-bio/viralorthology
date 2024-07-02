import os
from Bio import SeqIO
from tqdm import tqdm
from modules import commands
from modules.misc import get_n_genomes
from modules.misc import get_nseqs
from modules.misc import get_genomes_fasta
from modules.misc import get_file_list
from modules.blastp import make_db_protsingroup
from modules.misc import get_first_bioseq_fasta
from modules.misc import get_blastp_hits
from modules.misc import delete_tmp_files
from modules.misc import align_fasta_muscle
from modules.misc import delete_files
from modules.misc import combine_fastas
from modules.misc import get_ordered_files

def align_buildhmm():
    fastas = get_file_list()
    fastas = get_ordered_files(fastas)
    n_genomes = get_n_genomes()
    for fasta in fastas:
        n_genes = get_nseqs(fasta)
        if n_genomes == n_genes:
            continue
        align_fasta_muscle(fasta)
        output_a2m = fasta.replace('.fasta', '.a2m')
        commands.fasta_to_a2m(fasta.replace('.fasta', '.muscle'), output_a2m)
        commands.hmm_make(output_a2m, n_genomes)
        os.remove(output_a2m)

def get_hmmvshmm_score(g1, g2):
    report = f'{g1}-{g2}'
    commands.hmm_align(g1, g2, report)
    with open(report, encoding='utf-8') as report_op:
        score = 0
        for line in report_op:
            if line.strip().startswith('Probab='):
                score = float(line.strip().split()[0].split('=')[1])
                break
    os.remove(report)
    return score

def run_read_blastp(fasta):
    first_seq = get_first_bioseq_fasta(fasta)
    SeqIO.write(first_seq, 'query.fasta', 'fasta')
    commands.blastp('query.fasta', 'protDB_ingroup.fasta', '-evalue 10 -word_size 2')
    os.remove('query.fasta')
    hits = get_blastp_hits()
    groups_to_try = list(set([hit[1].replace('.fasta', '.hhm') for hit in hits])) #only the names, no evalues
    return groups_to_try

def get_compatible_groups():
    hmms = get_file_list('.hhm')
    final_groups = []
    already_checked = []
    print('Comparing HMMs...')
    for g1 in tqdm(hmms):
        if g1 in already_checked:
            continue
        already_checked.append(g1)
        compatible_groups = [g1]
        groups_to_try = run_read_blastp(g1)
        for g2 in groups_to_try:
            if g2 in already_checked:
                continue
            score = get_hmmvshmm_score(g1, g2)
            if score > 98: # MIN SCORE TO JOIN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                already_checked.append(g2)
                compatible_groups.append(g2)
        if len(compatible_groups) > 1:
            final_groups.append(compatible_groups)
    return final_groups

def combine_groups(compatible_groups):
    '''
    IN: list withs lists ([[fasta1, fasta2, fasta3], [fasta10, fasta12]])
    '''
    for group in compatible_groups:
        genomes = []
        group = [fasta.replace('.hhm', '.fasta') for fasta in group]
        for fasta in group:
            genomes += get_genomes_fasta(fasta)
        if len(genomes) == len(set(genomes)): #no repetead genomes, combine groups
            combine_fastas(group, group[0])
            delete_files([fasta.replace('.fasta', '.muscle') for fasta in group])

def hmmvshmm():
    align_buildhmm()
    make_db_protsingroup()
    compatible_groups = get_compatible_groups()
    combine_groups(compatible_groups)

def main():
    os.chdir('orthology_groups')

    hmmvshmm()
    delete_tmp_files(['.hhm'])
    os.remove('protDB_ingroup.fasta')

    os.chdir('..')
