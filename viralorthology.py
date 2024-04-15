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
from modules.prots_per_group import main as prots_per_group
from modules.merge_protDBs import main as merge_protDBs
from modules.synteny import main as synteny
from modules.check_filtered_genomes import main as check_filtered_genomes
from modules.update_db import main as update_db
from modules.download_seqs import main as download_seqs
from modules.check_dependencies import main as check_dependencies
from modules.kimura import main as kimura
from modules.misc import make_nseq_report
from modules.misc import delete_tmp_files_final

def split_params(usr_input):
    '''
    OUT: params[software] = [params]
    '''
    params = {}

    software = False
    for param in usr_input:
        if param in ['--proteinortho',
        '--orffinder',
        '--blastp',
        '--HMMsearch']:
            software = param.strip('-')
            params[software] = []
            continue
        if software:
            params[software].append(param)

    return params

def get_params(software, params):
    if software in params:
        soft_params = (' ').join(params[software])
        return soft_params
    return ''

def write_log(time, date, params):
    log = open('log.txt', 'w', encoding="utf-8")
    log.write(f'{date}\nElapsed time: {time}\n')
    log.write('Parameters:\n')
    for soft in params:
        params[soft] = (' ').join(params[soft])
        log.write(f'{soft}: {params[soft]}\n')
    log.close()

def check_params(proteinortho, blastp, orfinder, hmmsearch):
    '''
    checks if the user changed some prohibited parameter
    '''
    ok = True
    if proteinortho:
        if '--project' in proteinortho:
            ok = False
    if blastp:
        if '-query' in blastp or '-db' in blastp or '-out' in blastp:
            ok = False
    if orfinder:
        if '-in' in orfinder or '-out' in orfinder or '-outfmt' in orfinder:
            ok = False
    if not ok:
        print('Some prohibited parameter was changed')
        print('visit github for more information')
        sys.exit(1)

def main():
    start = time.time()

    usr_input = sys.argv[1:]
    if '-prots_per_group' in usr_input:
        prots_per_group()
        sys.exit(0)
    elif '-update_db' in usr_input:
        update_db()
        sys.exit(0)
    elif '-download_seqs' in usr_input:
        download_seqs()
        sys.exit(0)
    elif '-check_dependencies' in usr_input:
        check_dependencies()
        sys.exit(0)
    elif '-kimura' in usr_input:
        kimura()
        sys.exit(0)

    params_all = split_params(usr_input)
    params_po = get_params('proteinortho', params_all)
    params_blastp = get_params('blastp', params_all)
    params_orfinder = get_params('orffinder', params_all)
    params_HMMsearch = get_params('HMMsearch', params_all)
    check_params(params_po, params_blastp, params_orfinder, params_HMMsearch)

    #HERE COMES THE PIPELINE
    first_round = True
    while True:
        if first_round:
            split_data()
            filter_genomes()

        orfinder(params_orfinder)
        paralogs()
        make_protDB()

        if first_round:
            proteinortho(params_po)
            rename_groups()
            merge_groups()

            HMM_clean()
            HMMvsHMM()

            cleanDB()
            make_nseq_report('1-ProteinOrtho')

        HMM(params_HMMsearch)
        make_nseq_report('2-HMM')

        blastp('protDB.db', params_blastp)
        blastp('protDB_OF.db', params_blastp)
        make_nseq_report('3-Blastp')

        HMM(params_HMMsearch, True)
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
    write_log(deltatime, date, params_all)

if __name__ == '__main__':
    main()