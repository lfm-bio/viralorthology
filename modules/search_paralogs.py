import os
from Bio import SeqIO
from tqdm import tqdm
from modules import commands
from modules.misc import get_file_list
from modules.misc import get_blastp_hits
from modules.misc import get_bioseqs
from modules.misc import delete_tmp_files

def make_fastas_DB(seqs, proteome):
    '''
    makes blastdb with proteome and each query.fasta file
    '''
    commands.makeblastdb_prot(proteome)
    for seq in seqs:
        SeqIO.write(seq, f'{seq.id}.fasta', format='fasta-2line')
        fasta = f'{seq.id}.fasta'
        commands.makeblastdb_prot(fasta)

def search_paralogs(query, groups_paralogs, proteome):
    '''
    Searches for reciprocal paralogs of the query seq and adds them to the already found paralogs (groups_paralogs)
    '''
    paralogs = []
    commands.blastp(f'{query}.fasta', proteome, '-evalue 0.0000001')
    hits = get_blastp_hits()
    hits = [hit[1] for hit in hits] #only the IDs
    for hit in hits:
        commands.blastp(f'{hit}.fasta', f'{query}.fasta', '-evalue 0.0000001')
        hits = get_blastp_hits()
        if hits:
            paralogs.append(hit)
    paralogs.sort()

    if len(paralogs) > 1 and paralogs not in groups_paralogs:
        groups_paralogs.append(paralogs)
    return groups_paralogs

def keep_bigger_group(groups_paralogs):
    '''
    deletes groups that are inside other groups
    '''
    to_delete = []

    for n, g1 in enumerate(groups_paralogs):
        g1 = set(g1)
        for i, g2 in enumerate(groups_paralogs):
            if n == i:
                continue
            g2 = set(g2)

            if g1 == g2:
                to_delete.append(n) #group already is in the list

            if g1.issubset(g2):
                to_delete.append(n)

            if g2.issubset(g1):
                to_delete.append(i)

    if to_delete:
        for i in sorted(set(to_delete), reverse=True): #deletes from the end so indexes dont change
            del groups_paralogs[i]

    return groups_paralogs

def iteration(old):
    '''
    adds up groups that overlap
    '''
    new = []
    first = True

    while len(old) != len(new): #until nothing changes
        if not first:
            old = new
        first = False
        new = []
        for g1 in old:
            for g2 in old:
                if (set(g1)).intersection(set(g2)):
                    if sorted(set(g1+g2)) not in new:
                        new.append(sorted(set(g1 + g2)))
    return new

def clean_groups(groups_paralogs):
    final = []
    new_groups = iteration(groups_paralogs)

    #so it doesnt add reapeated groups
    for g in new_groups:
        if g not in final:
            final.append(g)

    for g in groups_paralogs:
        if g not in final:
            final.append(g)

    final = keep_bigger_group(final)
    return final

def delete_fastas(seqs_ids):
    for seq in seqs_ids:
        os.remove(f'{seq}.fasta')

def search_biggest_paralog(seqs, groups_paralogs):
    '''
    OUT: dict[biggest_prot] = [smaller_ones]
    '''
    paralogs_to_remove = {}
    for group in groups_paralogs:

        group = [seq for seq in seqs if seq.id in group]
        group.sort(key=lambda seq: len(seq.seq), reverse=True)
        group = [seq.id for seq in group]
        paralogs_to_remove[group[0]] = group[1:]

    return paralogs_to_remove

def make_new_files(proteome, prots_to_remove):
    '''
    removes small paralogs from proteome and makes paralogs/proteome-bigparalog.fasta with the smaller ones
    '''
    proteins_to_remove = []
    for big_paralog in prots_to_remove:
        small_paralogs = []
        for seq in SeqIO.parse(proteome, 'fasta'):
            if seq.id in prots_to_remove[big_paralog]: #if i want the biggest paralog to also be in the paralogs file add: or seq.id == big_paralog
                small_paralogs.append(seq)
                proteins_to_remove.append(seq.id)
        output_name = f'{proteome.replace(".fasta", "")}-{big_paralog}.fasta'
        SeqIO.write(small_paralogs, f'../paralogs/{output_name}', format='fasta-2line')

    #edits proteome
    to_new_proteome = []
    for seq in SeqIO.parse(proteome, 'fasta'):
        if seq.id not in proteins_to_remove:
            to_new_proteome.append(seq)

    SeqIO.write(to_new_proteome, 'newproteome.fasta', format='fasta-2line')
    os.remove(proteome)
    os.rename('newproteome.fasta', proteome)

def submain():
    proteomes = get_file_list()
    print('Searching for paralogs')
    for proteome in tqdm(proteomes):
        groups_paralogs = [] #paralogy groups that has been found [[g1], [g2],...]
        seqs = get_bioseqs(proteome) #seqs (biopython)
        make_fastas_DB(seqs, proteome) #makes blastdb with proteome and each protein
        seqs_ids = [seq.id for seq in seqs]

        for seq in seqs_ids:
            groups_paralogs = search_paralogs(seq, groups_paralogs, proteome)

        groups_paralogs = clean_groups(groups_paralogs)

        prots_to_remove = search_biggest_paralog(seqs, groups_paralogs) # dict[bigger] = [small1, small2, ...]
        make_new_files(proteome, prots_to_remove)

        delete_fastas(seqs_ids)

def main():
    if not os.path.isdir('paralogs'):
        os.mkdir('paralogs')
    os.chdir('proteomes')

    submain()
    delete_tmp_files(['.phr', '.pin', '.psq'])

    os.chdir('..')