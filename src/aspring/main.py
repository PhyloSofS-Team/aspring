import os
import glob
import subprocess # library to execute bash command line in python script
import argparse
from STEP4_dupraw2table import hhr2df
#from gene2ASRU.scripts.withoutDB.STEP4_dupraw2table import hhr2df
import sys
import pandas as pd

from aspring import __version__

def parse_args(args):
    parser = argparse.ArgumentParser(description='From Thoraxe outputs for a single query gene to its Alternative Splicing Repetitive Units')
    parser.add_argument('--gene',   dest='geneName', type=str, required=True, help='name of queried gene')
    parser.add_argument('--path_data',   type=str, required=True, help='path to dir containing Thoraxe outputs')
    parser.add_argument('--path_hhsuite',   type=str, required=True, help='path to dir of hhsuite ( YOURPATH/hhsuite or YOURPATH/hh-suite/build )')
    parser.add_argument('--len',   type=int, required=False, help='dont create profile for msa in which sequences are of length < X aa (def=5)',default=5)
    parser.add_argument('--id',     type=float, required=False, help='[0,100] maximum pairwise sequence identity (%%) (def=100)', default=100)
    parser.add_argument('--norealign',   type=int, required=True, help='bool, 1 if norealign else 0, do NOT realign displayed hits with Maximum Accuracy algorithm (MAC) (def=0)', default=0)
    parser.add_argument('--glo_loc',   type=int, required=True, help='bool, 1 if global else 0, use global/local alignment mode for searching/ranking (def=local)', default=0)
    parser.add_argument('--mact',   type=float, required=True, help= '[0,1[ posterior prob threshold for MAC realignment controlling greediness at alignment ends: 0:global >0.1:local (default=0.35)',default=0.35)
    parser.add_argument('--id_pair', type=float, required=True, help='[0,100] Identity percentage threshold between first sequence in msa of s-exon for each s-exon in a pair')
    parser.add_argument('--idCons_pair', type=float, required=True, help='[0,100] Identity percentage threshold between consensus sequence of msa of s-exon for each s-exon in a pair (should be equal to id_pair)')
    parser.add_argument('--pval', type=float, required=True, help='[0,1] p-value threshold for HMM-HMM alignment of a s-exons pair')
    parser.add_argument('--nbSpe', type=float, required=True, help='[1,10] minimum number of species in msa for s-exons in the pair')
    parser.add_argument('--cov', type=float, required=True, help='[0,1] Threshold for coverage of s-exon A and B in alignment of A and B')
    parser.add_argument('--version', action='version', version=f'aspring {__version__}')
    return parser.parse_args(args)

def run_pipeline(gene,path_data,path,msa_len,msa_id_threshold,re_align,glo_loc,mact,id_pair,idCons_pair,pval,nbSpe,cov):
    print("START STEP 1 : pre-processing : convert .fasta to .a2m")
    bashCommand=f"./ASPRIN/scripts/withoutDB/STEP1_preprocess_msas.py --gene {gene} --dataPATH {path_data} --hhsuitePATH {path} --len {msa_len}"
    subprocess.run([sys.executable]+bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 1")

    print("START STEP 2 : pre-processing : create HMM profiles")
    bashCommand=f"./ASPRIN/scripts/withoutDB/STEP2_hmm_maker.py --gene {gene} --path_data {path_data} --path {path} --id {msa_id_threshold}"
    subprocess.run([sys.executable]+bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 2")

    print("START STEP 3 : all-to-all pairwise HMM-HMM profile alignments")
    bashCommand=f"./ASPRIN/scripts/withoutDB/STEP3_hmm_aligner.py --gene {gene} --path_data {path_data} --path {path} --id {msa_id_threshold} --norealign {re_align} --glo_loc {glo_loc} --mact {mact}"
    subprocess.run([sys.executable]+bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 3")

    print("START STEP 4 : post-processing : parse alignment files and create corresponding table")
    os.makedirs(path_data+'/data', exist_ok=True) #make gene_dir in DupRaw, if folder already exist don't do anything
    os.makedirs(f"{path_data}/data/{gene}", exist_ok=True)
    #print(path_data+'/'+gene)
    if len(glob.glob(path_data+'/DupRaw/'+gene+'/*'))==0:
        sys.exit()
    if not os.path.exists(f'{path_data}/data/{gene}/{gene}_duplication_pairs.csv'):
        hhr2df(path_data+'/',gene).to_csv(f'{path_data}/data/{gene}/{gene}_duplication_pairs.csv',index=False)
    print("END STEP 4")

    print("START STEP 5 : post-processing : filter the table")
    bashCommand=f"./ASPRIN/scripts/withoutDB/STEP5_filter.py --gene {gene} --dataPATH {path_data} --id_pair {id_pair} --idCons_pair {idCons_pair} --pval {pval} --nbSpe {nbSpe} --cov {cov}"
    #if not os.path.exists(f"{path_data}/data/{gene}/{gene}_duplication_pairs_formated.csv"):
    subprocess.run([sys.executable]+bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 5")

    df=pd.read_csv(f'{path_data}/data/{gene}/{gene}_duplication_pairs_formated.csv')
    if df.shape[0]==0:
        print('There are no similar pairs of s-exons')
        sys.exit()

    print("START STEP 6 : post-processing : determine the type of similar s-exons in the context of AS")
    bashCommand=f"Rscript ./ASPRIN/scripts/withoutDB/STEP6_getStats4Massiv_gene.R --gene {gene} --path_data {path_data}"
    subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 6")

    #if not os.path.exists(path_data+"/data/"+gene+"/"+gene+"_duplication_pairs_filtered_valid_analEvents.csv"):
    #    sys.exit()
    print("START STEP 7 : post-processing : reformat table")
    bashCommand=f"./ASPRIN/scripts/withoutDB/STEP7_reformat.py --gene {gene} --dataPATH {path_data}"
    subprocess.run([sys.executable]+bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 7")

    print("START STEP 8 : create ASRUs and instances tables for query gene")
    bashCommand=f"./ASPRIN/scripts/withoutDB/STEP8_ASRU_2022_cleanedv3.py --gene {gene} --dataPATH {path_data}"
    subprocess.run([sys.executable]+bashCommand.split(), stdout=subprocess.PIPE)
    print("END STEP 8")

    print("START FINAL STEP : remove temporary files")
    bashCommand=f"./ASPRIN/scripts/withoutDB/STEP9_rm_tmpfiles.py --gene {gene} --path_data {path_data}"
    subprocess.run([sys.executable]+bashCommand.split(), stdout=subprocess.PIPE)
    print("END")

def run():
    args = parse_args(sys.argv[1:])

    gene = args.geneName
    path_data = args.path_data
    path = args.path_hhsuite
    msa_len = args.len
    msa_id_threshold = args.id # maximum pairwise sequence identity
    re_align = args.norealign
    glo_loc = args.glo_loc
    mact = args.mact
    id_pair = args.id_pair
    idCons_pair = args.idCons_pair
    pval = args.pval
    nbSpe = args.nbSpe
    cov = args.cov

    run_pipeline(gene,path_data,path,msa_len,msa_id_threshold,re_align,glo_loc,mact,id_pair,idCons_pair,pval,nbSpe,cov)


if __name__ == '__main__':
    run()