import pandas as pd
import glob
import os
import sys
import argparse
import warnings

from aspring import __version__


def read_fasta(fasta_file, keep_annotation=False):
    """Read sequences from fasta file.

    Parameters
    ----------
    fasta_file : str
        Name of fasta file to read.
    keep_annotation : boolean
        Determine is sequence id should contain annotation.

    Returns
    -------
    dict : dict[seq_id] -> seq
        Sequences indexed by sequence id.
    """

    if not os.path.exists(fasta_file):
        #raise InputFileError('Input file %s does not exist.' % fasta_file)
        print(f"{fasta_file} not found")

    if os.stat(fasta_file).st_size == 0:
        return {}

    try:

        if fasta_file.endswith('.gz'):
            file_f, file_mode = gzip.open, 'rt'
        else:
            file_f, file_mode = open, 'r'

        seqs = {}
        with file_f(fasta_file, file_mode) as f:

            for line in f.readlines():
                # skip blank lines
                if not line.strip():
                    continue

                if line[0] == '>':
                    if keep_annotation:
                        seq_id = line[1:-1]
                    else:
                        seq_id = line[1:].split(None, 1)[0]

                    seqs[seq_id] = []
                else:
                    seqs[seq_id].append(line.strip())

        for seq_id, seq in seqs.items():
            seqs[seq_id] = ''.join(seq).replace(' ', '')
    except:
        print(traceback.format_exc())
        print()
        print("[Error] Failed to process sequence file: " + fasta_file)
        sys.exit(1)

    return seqs


def hhr2df(path, gene):
    path_gene = path + '/DupRaw/' + gene
    hhr_files = glob.glob(path_gene + '/*.hhr')
    gene_df = []
    for aln in hhr_files:
        f = open(aln, 'r')
        hhr_info, gene_df_line = [], []
        for lines in f.readlines():
            if lines.strip() != '':
                hhr_info.append(lines.rstrip('\n').split())
        f.close()
        gene_df_line.append(hhr_info[0][1])  #s-exon query
        gene_df_line.append(hhr_info[8][1])  #s-exon hit
        gene_df_line.append(gene)  #gene
        gene_df_line.append(hhr_info[8][2])  #prob
        gene_df_line.append(hhr_info[8][3])  #e-val
        gene_df_line.append(hhr_info[8][4])  #p-val
        gene_df_line.append(hhr_info[8][5])  #score
        try:
            flag = False
            if '-' not in hhr_info[8][-3]:
                colQ = hhr_info[8][-2]
                flag = True
            else:
                colQ = hhr_info[8][-3]

            if flag == True and hhr_info[8][-1].rsplit('(')[0].split('-')[0] != '':
                colT = hhr_info[8][-1].rsplit('(')[0]
            else:
                colT = hhr_info[8][-2]
            taille_Q = float(hhr_info[1][-1])
            try:
                taille_T = float(hhr_info[8][-1])
            except ValueError:
                taille_T = float(hhr_info[16][-1].strip('()'))

            gene_df_line.append(colQ)
            gene_df_line.append(colT)
            gene_df_line.append(taille_Q)
            gene_df_line.append(taille_T)
            gene_df_line.append(hhr_info[11][4].rstrip('%').strip('Identities='))

            #compute id cons
            q = ""
            t = ""
            c = 0
            for l in hhr_info:
                if ' '.join(l).startswith("Q Consensus"):
                    q += l[3].upper()
                if ' '.join(l).startswith("T Consensus"):
                    t += l[3].upper()
            for k in range(min(len(q), len(t))):
                if q[k] == t[k]:
                    if q[k] != "~":
                        c = c + 1
            gene_df_line.append(str(100 * c / len(q)))

            gene_df_line.append(hhr_info[11][5].strip('Similarity='))
            gene_df_line.append(hhr_info[2][-1])  #nb species Q
            gene_df_line.append(
                len(
                    read_fasta(
                        f'{path}/{gene}/thoraxe/msa/msa_s_exon_{hhr_info[8][1]}.fasta'
                    )))
            gene_df.append(gene_df_line)
        # catch all the errors and show them as a warning
        except Exception as e:
            warnings.warn(f"hhr_info: {hhr_info}, error: {e}")

    df = pd.DataFrame(gene_df)
    df = df.rename(
        columns={
            0: 'S_exon_Q',
            1: 'S_exon_T',
            2: 'Gene',
            3: 'Prob',
            4: 'E-value',
            5: 'P-value',
            6: 'Score',
            7: 'Cols_Q',
            8: 'Cols_T',
            9: 'Length_Q',
            10: 'Length_T',
            11: 'Identities',
            12: 'IdCons',
            13: 'Similarity',
            14: 'NoSpecies_Q',
            15: 'NoSpecies_T'
        })

    return df


def parse_args(args):
    parser = get_arg_parser()
    return parser.parse_args(args)

def get_arg_parser():
    parser = argparse.ArgumentParser(
        description=
        'STEP 4 : Parses the alignment files and creates a table.'
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
    parser.add_argument(
        "--version",
        action="version",
        version=f"aspring {__version__}",
    )
    return parser


def dupraw2table(path_data, gene):
    os.makedirs(f'{path_data}/data', exist_ok=True)
    #make gene_dir in DupRaw, if folder already exist don't do anything
    os.makedirs(f"{path_data}/data/{gene}", exist_ok=True)
    if len(glob.glob(f'{path_data}/DupRaw/{gene}/*')) == 0:
        sys.exit()
    path_table = f"{path_data}/data/{gene}/{gene}_duplication_pairs.csv"
    if not os.path.exists(path_table):
        hhr2df(path_data, gene).to_csv(path_table, index=False)


def run():
    args = parse_args(sys.argv[1:])
    gene = args.geneName
    path_data = args.path_data
    dupraw2table(path_data, gene)


if __name__ == "__main__":
    run()
