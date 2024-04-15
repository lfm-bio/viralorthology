import os

def move_files():
    os.chdir('filtered')
    files = list(os.listdir(os.curdir))
    for xfile in files:
        os.replace(xfile, f'../{xfile}')
    os.chdir('..')
    os.rmdir('filtered')

def delete_move_files():
    for folder in ['genomes', 'orfeomes', 'proteomes']:
        os.chdir(folder)
        files = [xfile for xfile in os.listdir(os.curdir) if os.path.isfile(xfile)]
        for xfile in files:
            os.remove(xfile)
        move_files()
        os.chdir('..')

def main():
    delete_move_files()
    os.rename('n_seqs_groups.csv', 'first-round-n_seqs_groups.csv')
    os.rename('protDB_OF.db', 'first-round-protDB_OF.db')
    os.rename('protDB.db', 'first-round-protDB.db')