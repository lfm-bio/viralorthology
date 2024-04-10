import os

def check_filtered():
    os.chdir('filtered')

    files = list(os.listdir(os.curdir))
    return bool(files)

def main():
    os.chdir('genomes')
    filtered = check_filtered()
    os.chdir('../..')

    return filtered