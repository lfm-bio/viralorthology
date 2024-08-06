import os
from modules.misc import get_file_list
from modules.misc import get_nseqs
from modules.misc import get_all_genome_ids
from modules.misc import get_genomes_fasta
from modules.misc import get_ordered_files

def main():
    os.chdir('orthology_groups')

    report = open('../final_report.csv', 'w', encoding='utf-8')
    report.write('fasta,n_seqs,missing_genomes\n')
    genome_ids = get_all_genome_ids()

    fastas = get_file_list()
    fastas = get_ordered_files(fastas)
    for fasta in fastas:
        fasta_n_seqs = get_nseqs(fasta)
        fasta_genomes = get_genomes_fasta(fasta)
        missing_genomes = (',').join([genome for genome in genome_ids if genome not in fasta_genomes])
        report.write(f'{fasta},{fasta_n_seqs},{missing_genomes}\n')

    report.close()

    os.chdir('..')
