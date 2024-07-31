import os
from modules.misc import get_file_list
from modules.misc import get_ordered_files

def get_fasta_anotations(old_group):
    anotations = {}
    with open(old_group, encoding='utf-8') as fasta:
        for line in fasta:
            if line.startswith('>'):
                if '[protein=' in line:
                    anotation = line[line.find('[protein=')+9:line.find(']', line.find('[protein=')+9)]
                    anotation = anotation.replace(' ', '-').replace('/', '-').replace('(', '_').replace(')', '_').lower()
                    if 'unknown' in anotation:
                        continue
                    if anotation in anotations:
                        anotations[anotation] += 1
                    else:
                        anotations[anotation] = 1

    anotations = sorted([(n_seqs,an) for an,n_seqs in anotations.items()], reverse=True)
    anotations = [an for an in anotations if an[0] > 1] #more than one protein with that anotation

    if anotations:
        if anotations[0][1] == 'hypothetical-protein' and len(anotations) > 1:
            return anotations[1][1]
        return anotations[0][1]
    return 'hypothetical-protein'

def rename_files(new_names):
    for new_name, old_name in new_names.items():
        os.rename(old_name, new_name)

def rename_tmp_names():
    fastas = get_file_list()
    for n, fasta in enumerate(fastas):
        os.rename(fasta, f'tmp-{n}.fasta')

def rename_groups():
    rename_tmp_names()
    tmp_names = get_file_list()
    tmp_names = get_ordered_files(tmp_names)
    new_names = {}
    for old_file_name in tmp_names:
        new_name = get_fasta_anotations(old_file_name)
        if f'{new_name}.fasta' in new_names:
            n = 0
            while True:
                n += 1
                tmp_name = f'{new_name}-{n}'
                if f'{tmp_name}.fasta' not in new_names:
                    new_name = tmp_name
                    break
        new_name += '.fasta'
        new_names[new_name] = old_file_name
    rename_files(new_names)

def main():
    os.chdir('orthology_groups')

    rename_groups()

    os.chdir('..')
