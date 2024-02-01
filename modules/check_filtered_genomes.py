import os

def check_filtered():
    os.chdir('filtered')

    files = [xfile for xfile in os.listdir(os.curdir)]
    if not files:
        return False
    else:
        return True

def main():
    os.chdir('genomes')
    filtered = check_filtered()
    os.chdir('../..')

    if not filtered:
        return False
    else:
        return True