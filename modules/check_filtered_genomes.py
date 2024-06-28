import os

def main():

    def check_filtered():
        if os.path.isdir('filtered'):
            return True
        return False

    os.chdir('genomes')
    filtered = check_filtered()
    os.chdir('..')

    return filtered