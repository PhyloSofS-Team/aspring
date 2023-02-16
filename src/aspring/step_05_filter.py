import argparse
import pandas as pd
import sys

from aspring import __version__


def parse_args(args):
    parser = get_arg_parser()
    return parser.parse_args(args)

def get_arg_parser():
    parser = argparse.ArgumentParser(
        description='Filter the table to keep gene duplication pairs based on identity, coverage, p-value and number of species in the MSAs',
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
        '--id_pair',
        type=float,
        required=True,
        help=
        'Identity percentage threshold between first sequence in msa of s-exon for each s-exon in a pair'
    )
    parser.add_argument(
        '--idCons_pair',
        type=float,
        required=True,
        help=
        'Identity percentage threshold between consensus sequence of msa of s-exon for each s-exon in a pair'
    )
    parser.add_argument(
        '--pval',
        type=float,
        required=True,
        help='p-value threshold for HMM-HMM alignment of a s-exons pair')
    parser.add_argument(
        '--nbSpe',
        type=float,
        required=True,
        help='minimum number of species in msa for s-exons in the pair')
    parser.add_argument(
        '--cov',
        type=float,
        required=True,
        help='Threshold for coverage of s-exon A and B in alignment of A and B'
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"aspring {__version__}",
    )
    return parser


def filter(gene, path_data, id_pair, idCons_pair, pval, nbSpe, cov):

    intra_gene = pd.read_csv(path_data +
                             f'/data/{gene}/{gene}_duplication_pairs.csv')

    #filtrer la table
    intra_gene_filtered = intra_gene.loc[
        (intra_gene['Identities'] >= id_pair)
        & (intra_gene['IdCons'] >= idCons_pair) &
        (intra_gene['P-value'] <= pval) & (intra_gene['NoSpecies_Q'] >= nbSpe)
        & (intra_gene['NoSpecies_T'] >= nbSpe)]
    #!!!keep pairs where the alignement cover at least 80% of one of the s-exons
    ind_table = []
    for ind in intra_gene_filtered.index:
        row = intra_gene.iloc[ind]
        colsQ = int(row[7].split('-')[1]) - int(row[7].split('-')[0]) + 1
        coverageQ = colsQ / int(row[9])
        colsT = int(row[8].split('-')[1]) - int(row[8].split('-')[0]) + 1
        coverageT = colsT / int(row[10])
        if max(coverageQ, coverageT) >= 0.8:
            ind_table.append(ind)

    intra_gene = pd.DataFrame([intra_gene.iloc[ind] for ind in ind_table],
                              columns=list(intra_gene.columns))

    #reformattage pour appliquer script d'Elodie
    intra_gene['Cols'] = [
        int(i.split('-')[1]) - int(i.split('-')[0]) + 1
        for i in intra_gene['Cols_Q']
    ]
    intra_gene['Coverage_A'] = [
        100 * (int(row[7].split('-')[1]) - int(row[7].split('-')[0]) + 1) /
        int(row[9]) for row in intra_gene.to_numpy()
    ]
    intra_gene['Coverage_B'] = [
        100 * (int(row[8].split('-')[1]) - int(row[8].split('-')[0]) + 1) /
        int(row[10]) for row in intra_gene.to_numpy()
    ]
    intra_gene = intra_gene.drop(columns=[
        'Cols_Q', 'Cols_T', 'IdCons', 'NoSpecies_Q', 'NoSpecies_T', 'Score'
    ])
    intra_gene = intra_gene.drop(columns=["Identities", "Similarity"])
    intra_gene = intra_gene.rename(
        columns={
            "S_exon_Q": "S_exon_A",
            "S_exon_T": "S_exon_B",
            "E-value": "E_value",
            "P-value": "P_value",
            'Length_Q': 'Length_A',
            'Length_T': 'Length_B',
            'Identities': "Identity"
        })
    intra_gene = intra_gene[[
        'S_exon_A', 'S_exon_B', 'Prob', 'E_value', 'P_value', 'Cols', 'Gene',
        'Length_A', 'Length_B', 'Coverage_A', 'Coverage_B'
    ]]
    intra_gene.to_csv(
        f'{path_data}/data/{gene}/{gene}_duplication_pairs_formated.csv',
        index=0)

    fwrite = open(f"{path_data}/data/{gene}/{gene}_canonical_path.txt", 'w')
    f = open(path_data + "/" + gene + '/thoraxe/path_table.csv', 'r')
    tmp = [next(f) for x in range(2)]
    fwrite.write(gene + ' ' + tmp[1].split(',')[4] + '\n')
    fwrite.close()


def run():
    args = parse_args(sys.argv[1:])

    gene = args.geneName
    path_data = args.dataPATH
    id_pair = args.id_pair
    idCons_pair = args.idCons_pair
    pval = args.pval
    nbSpe = args.nbSpe
    cov = args.cov

    filter(gene, path_data, id_pair, idCons_pair, pval, nbSpe, cov)


if __name__ == "__main__":
    run()
