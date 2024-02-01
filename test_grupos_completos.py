import os
import sys

def search():
    fastas = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith('.fasta')]
    fastas.sort()
    try:
        n_seqs = int(sys.argv[1])
    except IndexError:
        print('\nUSO: grupos_completos.py n_seqs')
        quit()
    completos = 0
    for fasta in fastas:
        op = open(fasta)
        n = 0
        for line in op:
            if line.startswith('>'):
                n += 1
        op.close()

        if n == n_seqs:
            completos += 1
            print(fasta)
    print(f'Grupos completos: {completos}')

os.chdir('orthology_groups')
search()


