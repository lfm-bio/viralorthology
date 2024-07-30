import os
import copy
import statistics
from Bio import SeqIO
from modules.misc import get_file_list
from modules.misc import get_unique_gene_bioseq_by_id
from modules.misc import delete_tmp_files
from modules.misc import check_ortology_group
from modules.commands import blastp
from modules.commands import makeblastdb_prot
from modules.clean_protDB import cleanDB as clean_protDB

def check_synteny(conserved_windows, genes_in_order):

    def blastp_with_orthology_group(conserved_window, gene_found, synteny_window, ortho_groups):

        def blastp_hit():
            with open('blastp.results', encoding='utf-8') as results:
                for line in results:
                    if '***** No hits found *****' in line:
                        return False
            return True

        group_db = [group for group in conserved_window if group != gene_found][0] + '.fasta'
        makeblastdb_prot(group_db)
        querys = [gene for gene in synteny_window if gene not in ortho_groups]
        for query in querys:
            seq = get_unique_gene_bioseq_by_id(query)
            if seq is None:
                continue
            genome_id = seq.description.strip().split()[1]
            if check_ortology_group(genome_id, group_db) and seq:
                SeqIO.write(seq, 'query.fasta', 'fasta')
                blastp('query.fasta', group_db)
                if blastp_hit():
                    with open(group_db, 'a', encoding='utf-8') as fasta:
                        fasta.write(seq.format('fasta'))
                delete_tmp_files(['.pdb', '.phr', '.pin', '.pot', '.psq', '.ptf', '.pto'])
                os.remove('blastp.results')
                os.remove('query.fasta')
                clean_protDB()

    ortho_groups = get_file_list()
    ortho_groups = [fasta.replace('.fasta', '') for fasta in ortho_groups]
    for conserved_window in conserved_windows:
        for genome in genes_in_order:
            if not set(conserved_window).issubset(set(genes_in_order[genome])): #if that window is not found in that genome
                found = [gene for gene in conserved_window if gene in genes_in_order[genome]]
                if len(found) != 1: #if both or none were found i cant do anything
                    continue

                gene_found = found[0]
                found_pos_ingenome = genes_in_order[genome].index(gene_found)
                try:
                    synteny_window = genes_in_order[genome][found_pos_ingenome-1:found_pos_ingenome+2]
                except IndexError:
                    continue

                ok = False
                for gene in synteny_window:
                    if gene not in ortho_groups: #at least one has to be a unique gene
                        ok = True
                        break
                if ok:
                    blastp_with_orthology_group(conserved_window, gene_found, synteny_window, ortho_groups)

def get_conserved_windows(all_windows, n_genomes):
    output = open('../synteny.txt', 'w', encoding='utf-8')
    output.write('conservation synteny_window\n')
    conserved_windows = []
    for window in set(all_windows):
        perc = round(all_windows.count(window) / n_genomes * 100, 2)
        if 50 < perc < 100: #how many genome have that window? (if == 100% I wont find new genes)
            output.write(f'{perc}% {(" | ").join(window)}\n')
            conserved_windows.append(window)
    output.write('\n')
    output.close()

    return sorted(conserved_windows)

def get_genes(fasta, genes_no_order):
    '''
    OUT: dict[genome] = [(mid_pos, gene), (mid_pos, gene)...] NO ORDER
    '''

    def get_mid_pos(line):
        mid_pos = line[line.find('[location=')+10:line.find('[', line.find('[location=')+1)-2].strip('complement()')
        if line.startswith('>ORFINDER'):
            mid_pos = line[line.find('[location=')+10:line.find(']', line.find('[location=')+1)].strip('complement()')
        if mid_pos.startswith('join'):
            mid_pos = mid_pos[mid_pos.find('(')+1:].split(',')[0]
        mid_pos = mid_pos.replace('>', '').replace('<', '')
        mid_pos = [int(pos) for pos in mid_pos.split('..')]
        mid_pos = int(statistics.mean(mid_pos))
        return mid_pos

    no_order = copy.deepcopy(genes_no_order)
    gene = fasta.replace('.fasta', '')
    op_file = open(fasta, encoding='utf-8')
    for line in op_file:
        if line.startswith('>'):
            if fasta in ['../protDB.db', '../protDB_OF.db']:
                gene = line[1:].split()[0] # id instead of group name
            genome = line.split()[1]
            mid_pos = get_mid_pos(line)
            gene_data = (mid_pos, gene)
            if mid_pos:
                if genome in no_order:
                    no_order[genome].append(gene_data)
                else:
                    no_order[genome] = [gene_data]
    op_file.close()
    return no_order

def get_genes_in_order():

    def get_genes_no_order(fastas):
        '''
        OUT: dict[genome] = [(mid_pos, gene), (mid_pos, gene)...] NO ORDER
        '''
        genes_no_order = {}
        for fasta in fastas:
            genes_no_order = get_genes(fasta, genes_no_order)
        genes_no_order = get_genes('../protDB.db', genes_no_order)

        return genes_no_order

    def order_genes(genes):
        for genoma in genes:
            genes[genoma].sort()
        return genes

    def remove_midpos(genes_in_order):
        for genome in genes_in_order:
            genes_in_order[genome] = [gene[1] for gene in genes_in_order[genome]]

        return genes_in_order

    fastas = get_file_list()
    genes_no_order =  get_genes_no_order(fastas)
    genes_in_order = order_genes(genes_no_order)
    genes_in_order = remove_midpos(genes_in_order)

    return genes_in_order

def synteny():

    def get_windows(genes_in_order_onegenome):
        win_len = 2
        windows = []
        first, last = 0, win_len

        while True:
            win = genes_in_order_onegenome[first:last]
            if len(win) < win_len:
                break
            windows.append(tuple(win))
            first += 1
            last += 1

        return windows

    print('Synteny analysis...')
    genes_in_order = get_genes_in_order()
    n_genomes = len(genes_in_order)
    all_synteny_windows = []
    for ordered_genes in genes_in_order.values():
        all_synteny_windows += get_windows(ordered_genes)
    conserved_windows = get_conserved_windows(all_synteny_windows, n_genomes)
    check_synteny(conserved_windows, genes_in_order)

def main():
    os.chdir('orthology_groups')

    synteny()
    delete_tmp_files(['.pot', '.pin', '.psq', '.ptf', '.pto', '.phr', '.pdb'])

    os.chdir('..')
