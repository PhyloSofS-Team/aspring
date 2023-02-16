import os
import sys
import subprocess  # library to execute bash command line in python script
import glob
import argparse

from aspring import __version__


def parse_args(args):
    parser = get_arg_parser()
    return parser.parse_args(args)

def get_arg_parser():
    parser = argparse.ArgumentParser(
        description=
        'STEP 2 : Generates a Hidden Markov Model (HMM) profile for each s-exon.'
    )
    parser.add_argument('--gene',
                        dest='geneName',
                        type=str,
                        required=True,
                        help='name of gene')
    parser.add_argument(
        '--id',
        type=float,
        required=False,
        help='[0,100] maximum pairwise sequence identity (%%) (def=100)',
        default=100)
    parser.add_argument('--path_data',
                        type=str,
                        required=True,
                        help='path to dir containing Thoraxe outputs')
    #parser.add_argument('--M', type=float, required=False, help='[0,100] columns with fewer than X%% gaps are match states (def=50)', default=50)
    #parser.add_argument('--qsc', type=float, required=False, help='[0,100] minimum score per column with query (def=0)', default=0)
    #parser.add_argument('--len',   type=int, required=False, help='dont create profile for msa in which sequences are of length < X aa (def=5)', default=5)
    parser.add_argument(
        "--version",
        action="version",
        version=f"aspring {__version__}",
    )
    return parser


def hmm_maker(gene, msa_id_threshold, path_data):
    """Create HMM profiles for all s-exons of a chosen gene"""

    msa_folder = f'{path_data}/{gene}/thoraxe/msa'
    files = glob.glob(f"{msa_folder}/*")
    profiles = [profile for profile in files if '.hhm' in profile]
    a2ms = [msa for msa in files if '.a2m' in msa]

    if len(profiles) > 0:
        #not == len(msas) because cases where we won't create profiles (len(msa)<5aa)
        print('all HMM profiles are already created')
        exit()

    # Creating HMM profiles for gene

    for msa in a2ms:
        path = '/'.join(msa.split('/')[:-1])
        s_exon_name = msa.split(
            '/')[-1][:-4]  #-6 because it removes .a2m suffix
        s_exon_id = s_exon_name[11:]
        hmm_out = path + '/' + f'{s_exon_name}.hhm'
        if os.path.isfile(hmm_out):
            print('hmm_out is already created')
            continue
        bashCommand = f"hhmake -i {msa} -o {hmm_out} -v 0 -name {s_exon_id} -id {msa_id_threshold} -M a2m"
        subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)


def run():
    args = parse_args(sys.argv[1:])

    gene = args.geneName  # name of gene
    msa_id_threshold = args.id  # maximum pairwise sequence identity
    #msa_gap_threshold = args.M    # columns with fewer than X\% gaps are match states
    #msa_col_score = args.qsc    # columns with fewer than X\% gaps are match states
    path_data = args.path_data
    #msa_len = args.len

    hmm_maker(gene, msa_id_threshold, path_data)


if __name__ == '__main__':
    run()
