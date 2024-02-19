import os
from modules.misc import get_file_list
from modules.misc import get_ordered_files

def get_anotations(old_group):
    anotations = {}
    fasta = open(old_group)
    for line in fasta:
        if line.startswith('>'):
            if '[protein=' in line:
                anotation = line[line.find('[protein=')+9:line.find(']', line.find('[protein=')+9)]
                anotation = anotation.replace(' ', '-').replace('/', '-').replace('(', '_').replace(')', '_').lower()
                if anotation in anotations:
                    anotations[anotation] += 1
                else:
                    anotations[anotation] = 1
    fasta.close()
    anotations = sorted([(n_seqs,an) for an,n_seqs in anotations.items()], reverse=True)
    anotations = [an for an in anotations if an[0] > 1] #more than one protein with that anotation

    if anotations:
        if anotations[0][1] == 'hypothetical-protein' and len(anotations) > 1:
            return anotations[1][1]
        else:
            return anotations[0][1]
    else:
        return 'hypothetical-protein'

def rename_files(new_names):
    for new_name in new_names:
        os.rename(new_names[new_name], new_name)

def make_tmp_names():
    fastas = get_file_list()
    for n, fasta in enumerate(fastas):
        os.rename(fasta, f'tmp-{n}.fasta')

def rename_groups():
    make_tmp_names()
    old_names = get_file_list()
    old_names = get_ordered_files(old_names)
    new_names = {}
    for old_group in old_names:
        new_name = get_anotations(old_group)
        if f'{new_name}.fasta' in new_names:
            n = 0
            while True:
                n += 1
                tmp_name = f'{new_name}-{n}'
                if f'{tmp_name}.fasta' not in new_names:
                    new_name = tmp_name
                    break
        new_name += '.fasta'
        new_names[new_name] = old_group
    rename_files(new_names)

def main():
    os.chdir('orthology_groups')

    rename_groups()

    os.chdir('..')