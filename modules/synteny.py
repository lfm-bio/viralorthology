import os
import statistics
from modules.misc import get_file_list
import copy

def get_mid_pos(line):
    '''
    devuelve posicion central del gen
    '''
    mid_pos = line[line.find('[location=')+10:line.find('[', line.find('[location=')+1)-2].strip('complement()')
    if line.startswith('>ORFINDER'):
        mid_pos = line[line.find('[location=')+10:line.find(']', line.find('[location=')+1)].strip('complement()')
    if mid_pos.startswith('join'):
        mid_pos = mid_pos[mid_pos.find('(')+1:].split(',')[0]
    mid_pos = mid_pos.replace('>', '').replace('<', '')
    mid_pos = [int(pos) for pos in mid_pos.split('..')]
    mid_pos = int(statistics.mean(mid_pos))
    return mid_pos

def order_genes(genes):
    '''
    devuelve diccionario con lista de genes ordenados
    '''
    for genoma in genes:
        genes[genoma].sort()
    return genes

def get_genes(fasta, genes_no_order):
    '''
    OUT: dict[genome] = [(mid_pos, gene), (mid_pos, gene)...] NO ORDER
    '''
    no_order = copy.deepcopy(genes_no_order)
    gene = fasta.replace('.fasta', '')
    op_file = open(fasta)
    for line in op_file:
        if line.startswith('>'):
            if fasta == '../protDB.db' or fasta == '../protDB_OF.db':
                gene = line[1:].split()[0]
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

def get_genes_no_order(fastas):
    '''
    OUT: dict[genome] = [(mid_pos, gene), (mid_pos, gene)...] NO ORDER
    '''
    genes_no_order = {}
    for fasta in fastas:
        genes_no_order = get_genes(fasta, genes_no_order)
    genes_no_order = get_genes('../protDB.db', genes_no_order)
    genes_no_order_OF = get_genes('../protDB_OF.db', genes_no_order)

    return genes_no_order, genes_no_order_OF

def get_genes_in_order():
    '''
    devuelve dict genes[genoma] = [genes_ordenados] #(mid_pos, gene)
    '''
    fastas = get_file_list()
    genes_no_order, genes_no_order_OF =  get_genes_no_order(fastas)

    genes_in_order = order_genes(genes_no_order)
    genes_in_order_OF = order_genes(genes_no_order_OF)

    genes_in_order = remove_midpos(genes_in_order)
    genes_in_order_OF = remove_midpos(genes_in_order_OF)
    
    genes_in_order = genes_in_order_OF #NO OF GENES

    return genes_in_order, genes_in_order_OF

def remove_midpos(genes_in_order):

    for genome in genes_in_order:
        genes_no_midpos = []
        for gene in genes_in_order[genome]:
            genes_no_midpos.append(gene[1])
        genes_in_order[genome] = genes_no_midpos

    return genes_in_order

def get_windows(genes_in_order_onegenome, win_len):
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

def print_synteny(conserved_windows, genes_in_order):
    '''
    writes synteny.txt with the script output
    '''
    output = open('../synteny.txt', 'a')
    ortho_groups = get_file_list()
    ortho_groups = [fasta.replace('.fasta', '') for fasta in ortho_groups]
    for window in conserved_windows:
        for genome in genes_in_order:
            if not set(window).issubset(set(genes_in_order[genome])): #if that window is not found in that genome
                found = [gene for gene in window if gene in genes_in_order[genome]]
                if len(found) != 1: #if both or none were found i cant do anything
                    continue

                gene_found = found[0]
                found_pos_ingenome = genes_in_order[genome].index(gene_found)
                try:
                    todos = genes_in_order[genome][found_pos_ingenome-1:found_pos_ingenome+2]
                except:
                    continue

                ok = False
                for cada_gen in todos:
                    if cada_gen not in ortho_groups: #at least one has to be a unique gene
                        ok = True
                if ok:
                    output.write(f'conserved window: {(" | ").join(window)}\ngenome: {genome}\n{(" | ").join(todos)}\n\n')
    output.close()

def count_windows(windows_per_genome):
    output = open('../synteny.txt', 'w')
    output.write('conservation synteny_window\n')
    windows = []
    for genome in windows_per_genome:
        windows += windows_per_genome[genome]
    conserved_windows = []
    for window in set(windows):
        perc = round(windows.count(window) / len(windows_per_genome) * 100, 2)
        if perc > 40 and perc < 100: #how many genome have that window? (if 100% I wont find new genes)
            output.write(f'{perc}% {(" | ").join(window)}\n')
            conserved_windows.append(window)
    output.write('\n')
    output.close()
    return sorted(conserved_windows)

def get_synteny(genes_in_order, genes_in_order_OF):
    windows_per_genome = {genome : get_windows(genes_in_order[genome], 2) for genome in genes_in_order}
    conserved_windows = count_windows(windows_per_genome)
    print_synteny(conserved_windows, genes_in_order_OF)

def synteny():
    genes_in_order, genes_in_order_OF = get_genes_in_order() 
    get_synteny(genes_in_order, genes_in_order_OF)

def main():
    os.chdir('orthology_groups')

    synteny()

    os.chdir('..')
    print('check synteny.txt')