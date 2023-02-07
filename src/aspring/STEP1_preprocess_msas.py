import os
import glob
import subprocess # library to execute bash command line in python script
import argparse


#======================================== get custom user parameters ======================================

parser = argparse.ArgumentParser(description='Reformat fasta to a2m for all s-exons of a chosen gene. Necessary step before running hmm_maker !!')
parser.add_argument('--gene',   dest='geneName',    type=str, required=True, help='name of gene')
parser.add_argument('--dataPATH',   type=str, required=True, help='path to dir containing Thoraxe outputs')
parser.add_argument('--hhsuitePATH',   type=str, required=True, help='path to dir of hhsuite ( YOURPATH/hhsuite or YOURPATH/hh-suite/build )')
parser.add_argument('--len',   type=int, required=False, help='dont create profile for msa in which sequences are of length < X aa (def=5)',default=5)
args = parser.parse_args()

gene = args.geneName
path_data = args.dataPATH
path_hhsuite = args.hhsuitePATH
msa_len = args.len

#======================================== global variables ======================================

msa_folder='{}/{}/thoraxe/msa'.format(path_data,gene)
files=glob.glob("{}/*".format(msa_folder))
msas=[msa for msa in files if '.fasta' in msa]
a2ms=[file for file in files if '.a2m' in file]

#======================================== instructions ======================================

if len(a2ms) > 0: #not == len(msas) because cases where we won't create profiles (len(msa)<5aa)
    print('all MSAs are already converted')
    exit()

for msa in msas:
    f=open(msa,'r')
    beginning = [next(f).rstrip('\n') for x in range(2)] #get first sequence of msa (['>specie1', 'seq'])
    f.close()
    first_seq=beginning[1]
    if len(first_seq)<msa_len: #don't make profile for msa where len(seq) < 5 aa or user defined value
        continue
    else:
        s_exon_name=msa.split('/')[-1].rstrip('.fasta')
        bashCommand="{}/scripts/reformat.pl {} {}/{}/thoraxe/msa/{}.a2m".format(path_hhsuite, msa, path_data, gene, s_exon_name)
        subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)