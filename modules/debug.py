import os
from Bio import SeqIO

def write_n_seqs_per_group(stage):

    n_seqs_report = open('n_seqs_per_stage.txt', 'a', encoding='utf-8')
    os.chdir('orthology_groups')
    fastas = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith('.fasta')]
    fastas.sort()

    for fasta in fastas:
        if 'hypothetical-protein-' in fasta:
            continue
        n_seqs = 0
        for _ in SeqIO.parse(fasta, 'fasta'):
            n_seqs += 1
        n_seqs_report.write(f'{fasta},{stage},{n_seqs}\n')

    n_seqs_report.close()
    os.chdir('..')