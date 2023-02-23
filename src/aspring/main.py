import os
import subprocess  # library to execute bash command line in python script
import argparse
import sys
import pandas as pd

from aspring import __version__

def parse_args(args):
    parser = get_arg_parser()
    return parser.parse_args(args)

def get_arg_parser():
    parser = argparse.ArgumentParser(
        description=
        'From Thoraxe outputs for a single query gene to its Alternative Splicing Repetitive Units'
    )
    parser.add_argument('--gene',
                        dest='geneName',
                        type=str,
                        required=True,
                        help='name of queried gene')
    parser.add_argument('--path_data',
                        type=str,
                        required=False,
                        help='path to dir containing ThorAxe outputs (default: %(default)s)',
                        default=os.getcwd())
    parser.add_argument(
        '--path_hhsuite_scripts',
        type=str,
        required=False,
        help='path to the folder containing the scripts of hhsuite, by default it uses the value of the environment variable HHSUITE_SCRIPTS if it exists',
        default=os.environ.get('HHSUITE_SCRIPTS', ''))  # default value is empty string
    parser.add_argument(
        '--len',
        type=int,
        required=False,
        help=
        "don't create profiles for msas in which sequences are of length < len aa (default: %(default)s)",
        default=5)
    parser.add_argument(
        '--id',
        type=float,
        required=False,
        help='[0.0,100.0] maximum pairwise sequence identity (%%) (default: %(default)s)',
        default=100.0)
    parser.add_argument(
        '--norealign',
        type=int,
        required=False,
        help=
        'bool, 1 if norealign else 0, do NOT realign displayed hits with Maximum Accuracy algorithm (MAC) (default: %(default)s)',
        default=0)
    parser.add_argument(
        '--glo_loc',
        type=int,
        required=False,
        help=
        'bool, 1 if global else 0, use global/local alignment mode for searching/ranking (default: %(default)s)',
        default=0)
    parser.add_argument(
        '--mact',
        type=float,
        required=False,
        help=
        '[0.0,1.0] posterior prob threshold for MAC realignment controlling greediness at alignment ends: 0:global >0.1:local (default: %(default)s)',
        default=0.35)
    parser.add_argument(
        '--id_pair',
        type=float,
        required=False,
        help=
        '[0.0,100.0] Identity percentage threshold between first sequence in msa of s-exon for each s-exon in a pair (default: %(default)s)',
        default=50.0
    )
    parser.add_argument(
        '--idCons_pair',
        type=float,
        required=False,
        help=
        '[0.0,100.0] Identity percentage threshold between consensus sequence of msa of s-exon for each s-exon in a pair (should be equal to id_pair) (default: %(default)s)',
        default=50.0
    )
    parser.add_argument(
        '--pval',
        type=float,
        required=False,
        help='[0.0,1.0] p-value threshold for HMM-HMM alignment of a s-exons pair  (default: %(default)s)',
        default=0.001)
    parser.add_argument(
        '--nbSpe',
        type=int,
        required=False,
        help='[1,10] minimum number of species in msa for s-exons in the pair  (default: %(default)s)',
        default=2)
    parser.add_argument(
        '--cov',
        type=float,
        required=False,
        help=
        '[0.0,1.0] Threshold for coverage of s-exon A and B in alignment of A and B (default: %(default)s)',
        default=0.8
    )
    parser.add_argument('--version',
                        action='version',
                        version=f'aspring {__version__}')
    return parser


def check_if_executable_exists(executable):
    try:
        subprocess.run([executable, '--help'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        raise Exception(f'{executable} is not in the PATH. Please install it and try again.')


def run_pipeline(gene, path_data, path_hhsuite_scripts, msa_len,
                 msa_id_threshold, re_align, glo_loc, mact, id_pair,
                 idCons_pair, pval, nbSpe, cov):

    check_if_executable_exists('hhmake')
    check_if_executable_exists('hhalign')

    print("START STEP 1 : pre-processing : convert .fasta to .a2m")
    bashCommand = f"step_01_preprocess --gene {gene} --dataPATH {path_data} --path_hhsuite_scripts {path_hhsuite_scripts} --len {msa_len}"
    subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 1")

    print("START STEP 2 : pre-processing : create HMM profiles")
    bashCommand = f"step_02_hmm_maker --gene {gene} --path_data {path_data} --id {msa_id_threshold}"
    subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 2")

    print("START STEP 3 : all-to-all pairwise HMM-HMM profile alignments")
    bashCommand = f"step_03_hmm_aligner --gene {gene} --path_data {path_data} --id {msa_id_threshold} --norealign {re_align} --glo_loc {glo_loc} --mact {mact}"
    subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 3")

    print("START STEP 4 : parse alignments and create the corresponding table")
    bashCommand = f"step_04_gettable --gene {gene} --path_data {path_data}"
    subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 4")

    print("START STEP 5 : post-processing : filter the table")
    bashCommand = f"step_05_filter --gene {gene} --dataPATH {path_data} --id_pair {id_pair} --idCons_pair {idCons_pair} --pval {pval} --nbSpe {nbSpe} --cov {cov}"
    subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 5")

    df = pd.read_csv(
        os.path.join(path_data, 'data', gene, gene + '_duplication_pairs_formated.csv'))
    if df.shape[0] == 0:
        print('There are no similar pairs of s-exons')
        sys.exit()

    print(
        "START STEP 6 : post-processing : determine the type of similar s-exons in the context of AS"
    )
    bashCommand = f"step_06_stats --gene {gene} --path_data {path_data}"
    subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 6")

    print("START STEP 7 : post-processing : reformat table")
    bashCommand = f"step_07_reformat --gene {gene} --dataPATH {path_data}"
    subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 7")

    print("START STEP 8 : create ASRUs and instances tables for query gene")
    bashCommand = f"step_08_ASRUs --gene {gene} --dataPATH {path_data}"
    subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 8")

    print("START FINAL STEP : remove temporary files")
    bashCommand = f"step_09_clean --gene {gene} --path_data {path_data}"
    subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
    print("END")


def run():
    args = parse_args(sys.argv[1:])

    gene = args.geneName
    path_data = args.path_data
    path_hhsuite_scripts = args.path_hhsuite_scripts
    if not path_hhsuite_scripts:
        raise Exception('You need to specify the path to the hhsuite scripts using the --path_hhsuite_scripts argument or by setting the HHSUITE_SCRIPTS environment variable.')
    msa_len = args.len
    msa_id_threshold = args.id  # maximum pairwise sequence identity
    re_align = args.norealign
    glo_loc = args.glo_loc
    mact = args.mact
    id_pair = args.id_pair
    idCons_pair = args.idCons_pair
    pval = args.pval
    nbSpe = args.nbSpe
    cov = args.cov

    run_pipeline(gene, path_data, path_hhsuite_scripts, msa_len,
                 msa_id_threshold, re_align, glo_loc, mact, id_pair,
                 idCons_pair, pval, nbSpe, cov)


if __name__ == '__main__':
    run()