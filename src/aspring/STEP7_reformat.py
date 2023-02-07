import os
import glob
import argparse
import pandas as pd
import numpy as np

#======================================== get custom user parameters ======================================

parser = argparse.ArgumentParser(description='Reformat fasta to a2m for all s-exons of a chosen gene. Necessary step before running hmm_maker !!')
parser.add_argument('--gene',   dest='geneName',    type=str, required=True, help='name of gene')
parser.add_argument('--dataPATH',   type=str, required=True, help='path to dir containing Thoraxe outputs')
args = parser.parse_args()

gene = args.geneName
path_data = args.dataPATH

#======================================== instructions ======================================

eventsDup=pd.read_csv('{}/data/{}/{}_eventsDup.txt'.format(path_data,gene,gene))
tmp=pd.read_csv('{}/data/{}/{}_duplication_pairs.csv'.format(path_data,gene,gene)).drop_duplicates(['S_exon_Q','S_exon_T','Gene'])

print

new_cols=[]
for row_eve in eventsDup.to_numpy():
    vector=row_eve
    new_cols.append(tmp.loc[(tmp['S_exon_Q']==vector[1])&(tmp['S_exon_T']==vector[2])&(tmp['Gene']==vector[0])][['Cols_Q','Cols_T']].to_numpy())
    
new_cols=np.array(new_cols).reshape((-1,2))

eventsDup['ColA']=new_cols[:,0]
eventsDup['ColB']=new_cols[:,1]

eventsDup.to_csv('{}/data/{}/{}_eventsDup_withCols.txt'.format(path_data,gene,gene),index=0)