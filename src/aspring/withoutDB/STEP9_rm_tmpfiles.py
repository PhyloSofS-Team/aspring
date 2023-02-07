import os
import glob
import subprocess # library to execute bash command line in python script
import argparse




#======================================== get custom user parameters ======================================

parser = argparse.ArgumentParser(description='From Thoraxe outputs for a single query gene to its Alternative Splicing Repetitive Units')
parser.add_argument('--gene',   dest='geneName', type=str, required=True, help='name of queried gene')
parser.add_argument('--path_data',   type=str, required=True, help='path to dir containing Thoraxe outputs')
args = parser.parse_args()

gene = args.geneName
path_data = args.path_data

#======================================== instructions ======================================

os.remove(path_data+'/data/{}/{}_eventsDup.txt'.format(gene,gene))
os.remove(path_data+'/data/{}/{}_duplication_pairs_filtered_valid_analEvents.csv'.format(gene,gene))
os.remove(path_data+'/data/{}/{}_duplication_pairs_filtered_valid.csv'.format(gene,gene))
os.remove(path_data+'/data/{}/{}_duplication_pairs_formated.csv'.format(gene,gene))
os.remove(path_data+'/data/{}/{}_sexSizeEvents.csv'.format(gene,gene))
os.remove(path_data+'/data/{}/{}_canonical_path.txt'.format(gene,gene))

#remove .a2m files
for file in glob.glob(path_data+'/{}/thoraxe/msa/*.a2m'.format(gene)):
    os.remove(file)