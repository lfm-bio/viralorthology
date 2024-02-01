import os
import sys
from Bio import SeqIO


def get_all_genomes():
    genomes = []
    for seq in SeqIO.parse('genomes.fasta', 'fasta'):
        genomes.append(seq.description.split()[0])


    return genomes


fasta = sys.argv[1]

genomas_all = get_all_genomes()
os.chdir('orthology_groups')
genomas_fasta = []
fasta = open(fasta)
for line in fasta:
    if line.startswith('>'):
        genoma = line.strip().split()[1]
        genomas_fasta.append(genoma)
fasta.close()


for genoma in genomas_all:
    if genoma not in genomas_fasta:
        print(genoma)
