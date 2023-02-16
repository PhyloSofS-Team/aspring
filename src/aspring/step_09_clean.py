import os
import glob
import argparse
import sys

from aspring import __version__

def parse_args(args):
    parser = get_arg_parser()
    return parser.parse_args(args)

def get_arg_parser():
    parser = argparse.ArgumentParser(
        description=
        'Removes the intermediate files generated during the pipeline.'
    )
    parser.add_argument('--gene',
                        dest='geneName',
                        type=str,
                        required=True,
                        help='name of queried gene')
    parser.add_argument('--path_data',
                        type=str,
                        required=True,
                        help='path to dir containing Thoraxe outputs')
    parser.add_argument('--version',
                        action='version',
                        version=f'aspring {__version__}')
    return parser


def rm_tempfiles(gene, path_data):
    os.remove(f"{path_data}/data/{gene}/{gene}_eventsDup.txt")
    os.remove(
        f"{path_data}/data/{gene}/{gene}_duplication_pairs_filtered_valid_analEvents.csv"
    )
    os.remove(
        f"{path_data}/data/{gene}/{gene}_duplication_pairs_filtered_valid.csv")
    os.remove(f"{path_data}/data/{gene}/{gene}_duplication_pairs_formated.csv")
    os.remove(f"{path_data}/data/{gene}/{gene}_sexSizeEvents.csv")
    os.remove(f"{path_data}/data/{gene}/{gene}_canonical_path.txt")

    #remove .a2m files
    for file in glob.glob(path_data + '/{}/thoraxe/msa/*.a2m'.format(gene)):
        os.remove(file)


def run():
    args = parse_args(sys.argv[1:])
    gene = args.geneName
    path_data = args.path_data
    rm_tempfiles(gene, path_data)


if __name__ == '__main__':
    run()