#!/usr/bin/python3
import os
import sys
import time
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
from modules.check_filtered_genomes import main as check_filtered_genomes
from modules.download_seqs import main as download_seqs
from modules.check_dependencies import main as check_dependencies
from modules.kimura import main as kimura
from modules.misc import merge_protDBs
from modules.misc import make_nseq_report
from modules.misc import delete_tmp_files_final
from modules.misc import get_args
from modules.misc import delete_final_files
from modules.misc import write_log

def check_argv():
    if '-check_dependencies' in sys.argv:
        check_dependencies()
        sys.exit(0)
    if '-kimura' in sys.argv:
        kimura()
        sys.exit(0)
    if '-download_seqs' in sys.argv:
        download_seqs()
        sys.exit(0)

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

    check_argv()
    check_files()
    args = get_args()

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
            proteinortho(args.proteinortho)
            rename_groups()
            merge_groups()
            HMM_clean()
            HMMvsHMM()
            cleanDB()
            make_nseq_report('1-ProteinOrtho')

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

    merge_protDBs('protDB.db')
    merge_protDBs('protDB_OF.db')
    delete_tmp_files_final()
    delete_final_files()
    write_log(args, start)

if __name__ == '__main__':
    main()
