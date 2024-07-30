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
from modules.synteny import main as synteny
from modules.download_seqs import main as download_seqs
from modules.check_dependencies import main as check_dependencies
from modules.kimura import main as kimura
from modules.debug import write_n_seqs_per_group
from modules.misc import prepare_second_round
from modules.misc import check_filtered_genomes
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
            write_n_seqs_per_group('protein-ortho')
            merge_groups()
            HMM_clean()
            write_n_seqs_per_group('HMM-clean')
            HMMvsHMM()
            write_n_seqs_per_group('HMMvsHMM')
            clean_prot_db()

        HMM(args.HMMsearch)
        write_n_seqs_per_group('HMM-1')
        blastp('protDB.db', args.blastp)
        blastp('protDB_OF.db', args.blastp)
        write_n_seqs_per_group('BlastP')
        HMM(args.HMMsearch)
        write_n_seqs_per_group('HMM-2')

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
    write_n_seqs_per_group('Synteny')
    delete_tmp_files_final()
    delete_final_files()
    write_log(args, start)

if __name__ == '__main__':
    main()
