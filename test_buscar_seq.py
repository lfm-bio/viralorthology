#busca x_seq en todos los grupos de ortologia

import os
import sys

def search():
    fastas = [xfile for xfile in os.listdir(os.curdir) if xfile.endswith('.fasta')]
    fastas.sort()
    try:
        x_seq = sys.argv[1]
    except IndexError:
        print('\nUSO: buscar_seq.py ID')
        quit()

    for fasta in fastas:
        op = open(fasta)
        for line in op:
            if '[protein=' in line:
                anotation = line[line.find('[protein='):line.find(']', line.find('[protein='))]
                if x_seq.upper() in anotation.upper():
                    print(fasta)
                    break

        op.close()


    for fasta in fastas:
        op = open(fasta)
        for line in op:
            if line.startswith('>'):
                if x_seq.upper() in line.upper():
                    print(fasta)
                    break

        op.close()

os.chdir('orthology_groups')
search()


