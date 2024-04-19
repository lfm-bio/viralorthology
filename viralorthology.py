#!/usr/bin/python3
import os
import sys
import time
from datetime import timedelta
from datetime import datetime
from modules.split_data import main as split_data
from modules.orfinder import main as orfinder
from modules.proteinortho import main as proteinortho
from modules.make_protDB import main as make_protDB
from modules.search_paralogs import main as paralogs
from modules.blastp import main as blastp
from modules.HMM import main as HMM
from modules.clean_protDB import main as cleanDB
from modules.rename_groups import main as rename_groups
from modules.filter_genomes import main as filter_genomes
from modules.HMM_clean import main as HMM_clean
from modules.HMMVSHMM import main as HMMvsHMM
from modules.merge_groups import main as merge_groups
from modules.prepare_second_round import main as prepare_second_round
from modules.merge_protDBs import main as merge_protDBs
from modules.synteny import main as synteny
from modules.check_filtered_genomes import main as check_filtered_genomes
from modules.update_db import main as update_db
from modules.download_seqs import main as download_seqs
from modules.check_dependencies import main as check_dependencies
from modules.kimura import main as kimura
from modules.misc import make_nseq_report
from modules.misc import delete_tmp_files_final

class Args:
    def __init__(self):
        self.proteinortho = ''
        self.orffinder = ''
        self.blastp = ''
        self.HMMsearch = ''

    def get_all_params(self):
        return f'proteinortho: {self.proteinortho}\norfinder: {self.orffinder}\nblastp: {self.blastp}\nHMMsearch: {self.HMMsearch}'

    def __repr__(self):
        return f'{self.proteinortho} {self.orffinder} {self.blastp} {self.HMMsearch}'

def split_params(usr_input):
    '''
    OUT: params[software] = [params]
    '''
    args_dict = {}
    args = Args()

    software_list = [
        '--proteinortho',
        '--orffinder',
        '--blastp',
        '--HMMsearch'
        ]

    software = False
    for param in usr_input:
        if param in software_list:
            software = param.strip('-')
            args_dict[software] = []
            continue
        if software:
            args_dict[software].append(param)

    for software, params in args_dict.items():
        setattr(args, software.strip('-'), (' ').join(params))

    return args

def get_params(software, params):
    if software in params:
        soft_params = (' ').join(params[software])
        return soft_params
    return ''

def write_log(elapsed_time, date, args):
    with open('log.txt', 'w', encoding="utf-8") as log:
        log.write(f'{date}\nElapsed time: {elapsed_time}\n')
        log.write('Parameters:\n')
        log.write(args.get_all_params())

def check_params(args):
    '''
    checks if the user changed some prohibited parameter
    '''
    ok = True
    if args.proteinortho:
        if '--project' in args.proteinortho:
            ok = False
    if args.blastp:
        if '-query' in args.blastp or '-db' in args.blastp or '-out' in args.blastp:
            ok = False
    if args.orffinder:
        if '-in' in args.orffinder or '-out' in args.orffinder or '-outfmt' in args.orffinder:
            ok = False
    if not ok:
        print('Some prohibited parameter was changed')
        print('visit github for more information')
        sys.exit(1)

def check_argv(usr_input):
    if '-update_db' in usr_input:
        update_db()
        sys.exit(0)
    if '-download_seqs' in usr_input:
        download_seqs()
        sys.exit(0)
    if '-check_dependencies' in usr_input:
        check_dependencies()
        sys.exit(0)
    if '-kimura' in usr_input:
        kimura()
        sys.exit(0)

def only_on_first_round(params_po):
    proteinortho(params_po)
    rename_groups()
    merge_groups()

    HMM_clean()
    HMMvsHMM()

    cleanDB()
    make_nseq_report('1-ProteinOrtho')

def check_files():
    # ADD INFO
    if not os.path.isfile('genomes.fasta'):
        sys.exit(1)
    if not os.path.isfile('orfeomes.fasta'):
        sys.exit(1)
    if not os.path.isfile('proteomes.fasta'):
        sys.exit(1)

def main():
    start = time.time()

    usr_input = sys.argv[1:]
    check_argv(usr_input)

    check_files()

    args = split_params(usr_input)
    check_params(args)

    #HERE COMES THE PIPELINE
    first_round = True
    while True:
        if first_round:
            split_data()
            filter_genomes()

        orfinder(args.orffinder)
        paralogs()
        make_protDB()

        if first_round:
            only_on_first_round(args.proteinortho)

        HMM(args.HMMsearch)
        make_nseq_report('2-HMM')

        blastp('protDB.db', args.blastp)
        blastp('protDB_OF.db', args.blastp)
        make_nseq_report('3-Blastp')

        HMM(args.HMMsearch, True)
        make_nseq_report('4-HMM')

        rename_groups()
        merge_groups()

        if first_round:
            if not check_filtered_genomes(): #no filtered genomes
                break
            prepare_second_round()
            first_round = False
        else:
            break

    merge_protDBs()
    synteny()
    delete_tmp_files_final()

    os.system('rm -r genomes')
    os.system('rm -r orfeomes')
    os.system('rm -r proteomes')
    os.remove('n_seqs_groups.csv')
    if os.path.exists('first-round-n_seqs_groups.csv'):
        os.remove('first-round-n_seqs_groups.csv')
    os.remove('protDB_OF.db')

    end = time.time()
    deltatime = end - start
    deltatime = str(timedelta(seconds=deltatime))
    date = datetime.today().strftime('%Y-%m-%d')
    write_log(deltatime, date, args)

if __name__ == '__main__':
    main()