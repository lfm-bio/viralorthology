import os

def install_python_libs():
    os.system('pip install tqdm')
    os.system('pip install biopython')

def check():
    problems = []
    dependencies = {'muscle': 'muscle -h',
    'Entrez Direct': 'efetch -h',
    'BLAST': 'blastp -h',
    'ORFfinder': 'ORFfinder -h',
    'ProteinOrtho': 'proteinortho -h',
    'ProteinOrtho grab proteins': 'po_grab_proteins -h',
    'HMMER': 'hmmbuild -h',
    'hhsuite': 'hhmake -h',
    'clustalw': 'clustalw -help'}
    for k in dependencies:
        err = os.system(f'{dependencies[k]}')
        if err:
            if k == 'clustalw' and err == 256:
                continue
            problems.append(k)
    return problems

def main():
    install_python_libs()
    problems = check()
    print()
    if problems:
        print(f'Errors were found trying to run the following programs: {(" - ").join(problems)}')
    else:
        print('No problems were found.')