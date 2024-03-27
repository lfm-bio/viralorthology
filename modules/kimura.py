import os
from Bio import SeqIO
from modules.misc import get_file_list
from modules.misc import get_n_genomes
from modules.misc import get_nseqs
from modules.misc import get_bioseqs
from modules.misc import get_proteinid_genomeid
from modules.misc import fasta_to_fasta2line
from modules.misc import delete_files

def get_coregenes():
    '''
    OUT: list of orthology groups that have one protein for every genome
    '''
    fastas = get_file_list()
    n_genomes = get_n_genomes(True)
    coregenes = []
    for fasta in fastas:
        fasta_n_seqs = get_nseqs(fasta)
        if fasta_n_seqs == n_genomes:
            coregenes.append(fasta)
    return coregenes

def get_orfs():
    '''
    OUT: {genome_id: [orf1_bioseq, orf2_bioseq,...]}
    '''
    orfs = get_bioseqs('../orfeomes.fasta')
    orfs_per_genome = {}
    for seq in orfs:
        genome = seq.id[seq.id.find('|')+1:seq.id.find('.', seq.id.find('|'))+2]
        if genome not in orfs_per_genome:
            orfs_per_genome[genome] = [seq]
        else:
            orfs_per_genome[genome].append(seq)

    orfs_OF = get_bioseqs('../orfeomes_OF.fasta')
    for seq in orfs_OF:
        genome = seq.description.split()[1]
        if genome not in orfs_per_genome:
            orfs_per_genome[genome] = [seq]
        else:
            orfs_per_genome[genome].append(seq)

    return orfs_per_genome

def prot_to_orf(coregenes_fastas, orfs_per_genome):
    for fasta in coregenes_fastas:
        ids = get_proteinid_genomeid(fasta) # [(protID, genomeID), (protID, genomeID), ...]
        fasta_orfs = []
        for protID, genomeID in ids:
            for prot in orfs_per_genome[genomeID]:
                if protID in prot.id:
                    prot.id = genomeID
                    prot.description = genomeID
                    fasta_orfs.append(prot)
        fasta_orfs.sort(key=lambda x: x.id) #sort by genomeID
        SeqIO.write(fasta_orfs, f'../kimura/{fasta}', 'fasta')

def align_fastas():
    fastas = get_file_list()
    for fasta in fastas:
        os.system(f'clustalw -infile={fasta} -align -output=fasta')
        os.remove(fasta.replace('.fasta', '.dnd'))
        fasta_to_fasta2line(fasta)

def join_fastas():
    megaprots = {}
    fastas = get_file_list()
    for fasta in fastas:
        for seq in SeqIO.parse(fasta, 'fasta'):
            if seq.id not in megaprots:
                megaprots[seq.id] = str(seq.seq)
            else:
                megaprots[seq.id] += str(seq.seq)
    delete_files(fastas)
    finalfasta = open('final.fasta', 'w', encoding='utf-8')
    for genomeID in megaprots:
        finalfasta.write(f'>{genomeID}\n{megaprots[genomeID]}\n')
    finalfasta.close()

def main():
    os.mkdir('kimura')
    os.chdir('orthology_groups')

    orfs_per_genome = get_orfs()
    coregenes_fastas = get_coregenes()
    prot_to_orf(coregenes_fastas, orfs_per_genome)
    os.chdir('../kimura')
    align_fastas()
    join_fastas()

    os.chdir('..')