import sys
import argparse
import pandas as pd
import numpy as np

from aspring import __version__

def parse_args(args):
    parser = get_arg_parser()
    return parser.parse_args(args)

def get_arg_parser():
    parser = argparse.ArgumentParser(
        description=
        'Reformat the previous outputs to add the information about the duplicated regions.'
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
    parser.add_argument('--version',
                        action='version',
                        version=f'aspring {__version__}')
    return parser


def reformate(gene, path_data):
    eventsDup = pd.read_csv(f'{path_data}/data/{gene}/{gene}_eventsDup.txt')
    tmp = pd.read_csv(f'{path_data}/data/{gene}/{gene}_duplication_pairs.csv'
                      ).drop_duplicates(['S_exon_Q', 'S_exon_T', 'Gene'])

    print

    new_cols = []
    for row_eve in eventsDup.to_numpy():
        vector = row_eve
        new_cols.append(
            tmp.loc[(tmp['S_exon_Q'] == vector[1])
                    & (tmp['S_exon_T'] == vector[2])
                    & (tmp['Gene'] == vector[0])][['Cols_Q',
                                                   'Cols_T']].to_numpy())

    new_cols = np.array(new_cols).reshape((-1, 2))

    eventsDup['ColA'] = new_cols[:, 0]
    eventsDup['ColB'] = new_cols[:, 1]

    eventsDup.to_csv(f'{path_data}/data/{gene}/{gene}_eventsDup_withCols.txt',
                     index=0)


def run():
    args = parse_args(sys.argv[1:])
    gene = args.geneName
    path_data = args.dataPATH
    reformate(gene, path_data)


if __name__ == '__main__':
    run()