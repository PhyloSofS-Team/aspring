import os
import subprocess # library to execute bash command line in python script
import glob
import argparse

#======================================== get custom user parameters ======================================

parser = argparse.ArgumentParser(description='From pre-proccessed ThorAxe outputs, create HMM profiles for all s-exons of a chosen gene. Pre-processing step of converting .fasta to .a2m is mandatory !')
parser.add_argument('--gene',   dest='geneName',    type=str, required=True, help='name of gene')
parser.add_argument('--id',     type=float, required=False, help='[0,100] maximum pairwise sequence identity (%%) (def=100)', default=100)
parser.add_argument('--path_data',   type=str, required=True, help='path to dir containing Thoraxe outputs')
#parser.add_argument('--M', type=float, required=False, help='[0,100] columns with fewer than X%% gaps are match states (def=50)', default=50)
#parser.add_argument('--qsc', type=float, required=False, help='[0,100] minimum score per column with query (def=0)', default=0)
#parser.add_argument('--len',   type=int, required=False, help='dont create profile for msa in which sequences are of length < X aa (def=5)', default=5)
parser.add_argument('--path',   type=str, required=True, help='path to bin dir of hhsuite ( YOURPATH/hhsuite or YOURPATH/hh-suite/build )')  #change it to hh-suite for users that didn't downloaded it from https://mmseqs.com/hhsuite/
args = parser.parse_args()

gene = args.geneName     # name of gene
msa_id_threshold = args.id # maximum pairwise sequence identity
#msa_gap_threshold = args.M    # columns with fewer than X\% gaps are match states
#msa_col_score = args.qsc    # columns with fewer than X\% gaps are match states
path_data = args.path_data
hhsuitebin_path = args.path
#msa_len = args.len


#======================================== global variables ======================================

msa_folder='{}/{}/thoraxe/msa'.format(path_data,gene)
files=glob.glob("{}/*".format(msa_folder))
profiles=[profile for profile in files if '.hhm' in profile]
a2ms=[msa for msa in files if '.a2m' in msa]

#======================================== instructions ======================================

if len(profiles) > 0: #not == len(msas) because cases where we won't create profiles (len(msa)<5aa)
	print('all HMM profiles are already created')
	exit()

# Creating HMM profiles for gene

for msa in a2ms:
	path='/'.join(msa.split('/')[:-1])
	s_exon_name=msa.split('/')[-1][:-4] #-6 because it removes .a2m suffix
	s_exon_id=s_exon_name[11:]
	hmm_out=path+'/'+'{}.hhm'.format(s_exon_name)
	if os.path.isfile(hmm_out):
		print('hmm_out is already created')
		continue
	bashCommand = "{}/bin/hhmake -i {} -o {} -v 0 -name {} -id {} -M a2m".format(hhsuitebin_path, msa, hmm_out, s_exon_id, msa_id_threshold)
	subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
