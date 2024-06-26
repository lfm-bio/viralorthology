#!/usr/bin/python3
import sys
import time
from modules.split_data import main as split_fastas
from modules.orfinder import main as orfinder
from modules.proteinortho import main as proteinortho
from modules.make_protDB import main as make_protdb
from modules.search_paralogs import main as search_paralogs
from modules.blastp import main as blastp
from modules.HMM import main as HMM
from modules.clean_protDB import main as clean_prot_db
from modules.rename_groups import main as rename_groups
from modules.filter_genomes import main as filter_similar_genomes
from modules.HMM_clean import main as HMM_clean
from modules.HMMVSHMM import main as HMMvsHMM
from modules.merge_groups import main as merge_groups
from modules.prepare_second_round import main as prepare_second_round
from modules.synteny import main as synteny
from modules.check_filtered_genomes import main as check_filtered_genomes
from modules.download_seqs import main as download_seqs
from modules.check_dependencies import main as check_dependencies
from modules.kimura import main as kimura
from modules.misc import merge_protDBs
from modules.misc import delete_tmp_files_final
from modules.misc import get_args
from modules.misc import delete_final_files
from modules.misc import write_log
from modules.misc import check_files

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

def main():
    start = time.time()

    check_argv()
    check_files()
    args = get_args()

    #HERE COMES THE PIPELINE
    first_round = True
    while True:
        if first_round:
            split_fastas()
            filter_similar_genomes()

        orfinder(args.orffinder)
        search_paralogs()
        make_protdb()

        if first_round:
            proteinortho(args.proteinortho)
            rename_groups()
            merge_groups()
            HMM_clean()
            HMMvsHMM()
            clean_prot_db()

        HMM(args.HMMsearch)
        blastp('protDB.db', args.blastp)
        blastp('protDB_OF.db', args.blastp)
        HMM(args.HMMsearch)

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
    synteny()
    delete_tmp_files_final()
    delete_final_files()
    write_log(args, start)

if __name__ == '__main__':
    main()
