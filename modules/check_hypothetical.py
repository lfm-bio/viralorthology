import os
from modules.misc import get_file_list

def get_anotations(fasta):
    anotations = []
    fasta_op = open(fasta)
    for line in fasta_op:
        if line.startswith('>'):
            anotation = line[line.find('[protein=')+9:line.find(']', line.find('[protein=')+9)]
            anotation = anotation.replace(' ', '-').replace('/', '-').replace('(', '_').replace(')', '_').lower()
            if anotation != 'hypothetical-protein' and anotation != 'unknown' and not anotation.startswith('orf'):
                anotations.append(anotation)
    fasta_op.close()
    return anotations

def check():
    output = open('../check_hypothetical-proteins.txt', 'w')
    output.write(f'{"ortology_group":30s} hypothetical_group\n')
    fastas = get_file_list()
    hyphos = [fasta for fasta in fastas if fasta.startswith('hypothetical-protein')]
    no_hyphos = [fasta.replace('.fasta', '') for fasta in fastas if not fasta.startswith('hypothetical-protein')]
    for hypho in hyphos:
        anotations = get_anotations(hypho)
        for anotation in set(anotations):
            if anotation in no_hyphos:
                output.write(f'{anotation:30s} {hypho}\n')
    output.close()

def main():
    os.chdir('orthology_groups')

    check()
    print('Check check_hypothetical-proteins.txt')

    os.chdir('..')

if __name__ == '__main__':
    main()