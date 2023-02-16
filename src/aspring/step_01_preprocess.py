import os
import sys
import glob
import subprocess  # library to execute bash command line in python script
import argparse

from aspring import __version__

def parse_args(args):
    parser = get_arg_parser()
    return parser.parse_args(args)

def get_arg_parser():
    parser = argparse.ArgumentParser(
        description=
        'STEP 1 : Reformat s-exons fasta files to a2m'
    )
    parser.add_argument('--gene',
                        dest='geneName',
                        type=str,
                        required=True,
                        help='name of gene')
    parser.add_argument('--dataPATH',
                        type=str,
                        required=True,
                        help='path to dir containing Thoraxe outputs')
    parser.add_argument(
        '--path_hhsuite_scripts',
        type=str,
        required=True,
        help='path to the folder containing the scripts of hhsuite')
    parser.add_argument(
        '--len',
        type=int,
        required=False,
        help=
        'dont create profile for msa in which sequences are of length < X aa (def=5)',
        default=5)
    parser.add_argument(
        "--version",
        action="version",
        version=f"aspring {__version__}",
    )
    return parser


def preprocess_msas(gene, path_data, path_hhsuite, msa_len):
    msa_folder = os.path.join(path_data, gene, 'thoraxe', 'msa')
    files = glob.glob(f"{msa_folder}/*")
    msas = [msa for msa in files if '.fasta' in msa]
    a2ms = [file for file in files if '.a2m' in file]

    if len(a2ms) > 0:
        #not == len(msas) because cases where we won't create profiles (len(msa)<5aa)
        print('all MSAs are already converted')
        exit()

    for msa in msas:
        f = open(msa, 'r')
        beginning = [next(f).rstrip('\n') for x in range(2)]
        #get first sequence of msa (['>specie1', 'seq'])
        f.close()
        first_seq = beginning[1]
        if len(first_seq) < msa_len:
            #don't make profile for msa where len(seq) < 5 aa or user defined value
            continue
        else:
            s_exon_name = msa.split('/')[-1].rstrip('.fasta')
            a2m_file = os.path.join(path_data, gene, 'thoraxe', 'msa',
                                    f"{s_exon_name}.a2m")
            bashCommand = f"{os.path.join(path_hhsuite, 'reformat.pl')} {msa} {a2m_file}"
            subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)


def run():
    """Entry point for the application script"""
    args = parse_args(sys.argv[1:])

    gene = args.geneName
    path_data = args.dataPATH
    path_hhsuite = args.path_hhsuite_scripts
    msa_len = args.len

    preprocess_msas(gene, path_data, path_hhsuite, msa_len)


if __name__ == "__main__":
    run()
