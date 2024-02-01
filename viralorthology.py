#!/usr/bin/python3
import modules.split_data as split_data
import modules.orfinder as orfinder
import modules.proteinortho as proteinortho
import modules.make_protDB as make_protDB
import modules.search_paralogs as paralogs
import modules.blastp as blastp
import modules.HMM as HMM
import modules.clean_protDB as cleanDB
import modules.rename_groups as rename_groups
import modules.filter_genomes as filter_genomes
import modules.HMM_clean as HMM_clean
import modules.HMMVSHMM as HMMvsHMM
import modules.merge_groups as merge_groups
import modules.prepare_second_round as prepare_second_round
import modules.prots_per_group as prots_per_group
import modules.merge_protDBs as merge_protDBs
import modules.try_to_merge as try_to_merge
import modules.synteny as synteny
import modules.check_filtered_genomes as check_filtered_genomes
import modules.update_db as update_db
import modules.download_seqs as download_seqs
import modules.check_dependencies as check_dependencies
import modules.kimura as kimura
from modules.misc import make_nseq_report
from datetime import timedelta
from datetime import datetime
import os
import sys
import time

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
    else:
        return False

def write_log(time, date, params):
    log = open('log.txt', 'w')
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
        quit()

def main():
    start = time.time()

    usr_input = sys.argv[1:]
    if '-prots_per_group' in usr_input:
        prots_per_group.main()
        quit()
    elif '-try_to_merge_groups' in usr_input:
        fasta_list = usr_input[1:]
        try_to_merge.main(fasta_list)
        quit()
    elif '-update_db' in usr_input:
        update_db.main()
        quit()
    elif '-download_seqs' in usr_input:
        download_seqs.main()
        quit()
    elif '-check_dependencies' in usr_input:
        check_dependencies.main()
        quit()
    elif '-kimura' in usr_input:
        kimura.main()
        quit()

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
            split_data.main()
            filter_genomes.main()

        orfinder.main(params_orfinder)
        paralogs.main()
        make_protDB.main()

        if first_round:
            proteinortho.main(params_po)
            rename_groups.main()
            merge_groups.main()

            HMM_clean.main()
            HMMvsHMM.main()

            cleanDB.main()
            make_nseq_report('1-ProteinOrtho')

        HMM.main(params_HMMsearch)
        make_nseq_report('2-HMM')

        blastp.main('protDB.db', params_blastp)
        blastp.main('protDB_OF.db', params_blastp)
        make_nseq_report('3-Blastp')

        HMM.main(params_HMMsearch, True)
        make_nseq_report('4-HMM')

        rename_groups.main()
        merge_groups.main()

        if first_round:
            if not check_filtered_genomes.main(): #no filtered genomes
                break
            prepare_second_round.main()
            first_round = False
        else:
            break

    merge_protDBs.main()
    synteny.main()

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