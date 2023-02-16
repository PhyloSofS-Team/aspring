import os
import sys
import subprocess  # library to execute bash command line in python script
import glob
import itertools
import argparse
import datetime

from aspring import __version__


def parse_args(args):
    parser = get_arg_parser()
    return parser.parse_args(args)

def get_arg_parser():
    parser = argparse.ArgumentParser(
        description=
        'STEP 3 : HMM-HMM alignment of all the s-exons combinations.'
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
        help=
        '[0,100] maximum pairwise sequence identity (%%) applied to query MSA, template MSA, and result MSA (def=100)',
        default=100)
    parser.add_argument('--path_data',
                        type=str,
                        required=True,
                        help='path to dir containing Thoraxe outputs')
    parser.add_argument(
        '--norealign',
        type=int,
        required=True,
        help=
        'bool, 1 if norealign else 0, do NOT realign displayed hits with Maximum Accuracy algorithm (MAC) (def=0)',
        default=0)
    parser.add_argument(
        '--glo_loc',
        type=int,
        required=True,
        help=
        'bool, 1 if global else 0, use global/local alignment mode for searching/ranking (def=local)',
        default=0)
    parser.add_argument(
        '--mact',
        type=float,
        required=True,
        help=
        '[0,1[ posterior prob threshold for MAC realignment controlling greediness at alignment ends: 0:global >0.1:local (default=0.35)',
        default=0.35)
    parser.add_argument(
        "--version",
        action="version",
        version=f"aspring {__version__}",
    )
    return parser


def hmm_aligner(gene, msa_id_threshold, path_data, re_align,
                glo_loc, mact):
    """
	HMM-HMM profile alignment of all s-exons combinations for a chosen gene
	"""

    print('start', datetime.datetime.now().time())

    dico_gloloc = {0: '-local', 1: '-global'}
    dico_align = {0: '', 1: '-norealign'}

    msa_folder = f'{path_data}/{gene}/thoraxe/msa'
    profiles = glob.glob(f"{msa_folder}/*.hhm")
    #profiles=[profile for profile in files if '.hhm' in profile]

    if len(profiles) == 0:
        print('You need to generate the HMM profiles !!')
        exit()

    aln_folder = f'{path_data}/DupRaw'  #where the alignments will be located
    gene_aln_folder = aln_folder + '/' + gene

    os.makedirs(
        gene_aln_folder, exist_ok=True
    )  #make gene_dir in DupRaw, if folder already exist don't do anything

    for (hmm_a, hmm_b) in itertools.combinations(profiles, 2):
        name_a = hmm_a.split('/')[-1].rstrip('.hhm')
        name_b = hmm_b.split('/')[-1].rstrip('.hhm')
        id_a, id_b = name_a[11:], name_b[11:]  #name is msa_s_exon_X
        hhr_out = gene_aln_folder + '/' + id_a + '.' + id_b + '.hhr'
        bashCommand = f'hhalign -i {hmm_a} -t {hmm_b} -o {hhr_out} -v 0 -id {msa_id_threshold} {dico_align[re_align]} {dico_gloloc[glo_loc]} -mact {mact}'
        if not os.path.exists(hhr_out):
            subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)

    print(gene, 'end', datetime.datetime.now().time())


def run():
    args = parse_args(sys.argv[1:])

    gene = args.geneName
    msa_id_threshold = args.id
    path_data = args.path_data
    re_align = args.norealign
    glo_loc = args.glo_loc
    mact = args.mact

    hmm_aligner(gene, msa_id_threshold, path_data, re_align,
                glo_loc, mact)


if __name__ == '__main__':
    run()
