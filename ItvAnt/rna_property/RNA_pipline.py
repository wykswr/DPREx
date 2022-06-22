#!/home/fancong/anaconda3/bin/python3.8
import sys
import pandas as pd
import numpy as np
import os
import json
from subprocess import Popen, PIPE
#for fan calculation
from pathlib import Path
import annot_Rloop
import annot_BER
import get_bed
import score_RNA_pair

#calculation paths
with open(Path(__file__).parent / 'fan_path.json') as handle:
    config = json.load(handle)



#read input
fn_input=sys.argv[1]
pd_loc=pd.read_csv(fn_input,header=None,names=['chrom','start','end'],sep='\t')
pd_loc.iloc[:,[1,2]]=pd_loc.iloc[:,[1,2]].astype('int64')
#calculateRloop
l_rloop=pd_loc.index.map(lambda x:annot_Rloop.get_olp2(x,pd_loc))
#print(l_rloop)

#calculate BER, DNA_repair factors
df_BER_res=annot_BER.cal_sum1(pd_loc)
#print(df_BER_res)

#transfer to sequence
df1,df2=get_bed.new_pos1(pd_loc)
dir_out=config['pair_dir']
df1.to_csv(dir_out+'new_bed/df1.bed',index=False,header=False,sep='\t')
df2.to_csv(dir_out+'new_bed/df2.bed',index=False,header=False,sep='\t')
bedtools_run=config['bedtools']
cmd=bedtools_run+' getfasta -fi '+config['seq_for_bedtools']+' -bed '+dir_out+'new_bed/df1.bed -fo '+dir_out+'get_seq/df1.fa'
os.system(cmd)
cmd=bedtools_run+' getfasta -fi '+config['seq_for_bedtools']+' -bed '+dir_out+'new_bed/df2.bed -fo '+dir_out+'get_seq/df2.fa'
os.system(cmd)
#cal RNA_pair, fuck that! use the shell!
os.system('./run_GRASP.sh')
fn_df1=dir_out+'gf_grass/df1_upper.out'
fn_df2=dir_out+'gf_grass/df2_upper.out'
res_pair1=score_RNA_pair.read_score(fn_df1,0) #start
res_pair2=score_RNA_pair.read_score(fn_df2,1) #end
#print(res_pair1)
#cal free energy
os.system('./RNA_FE.sh {}get_seq/df1.fa'.format(dir_out))
os.system('./RNA_FE.sh {}get_seq/df2.fa'.format(dir_out))
dir_out=config['FE_dir']
res_FE1=pd.read_csv('{}df1.energ'.format(dir_out),sep='\t|-|_',header=None,names=['chrom','start','energ'],usecols=[0,2,4],engine='python')
res_FE2=pd.read_csv('{}df2.energ'.format(dir_out),sep='\t|-|_',header=None,names=['chrom','end','energ'],usecols=[0,1,4],engine='python')
# print(res_FE1)

#combine result
df_RNA_res=pd_loc.loc[:,:]
df_RNA_res['Rloop']=l_rloop
df_RNA_res=pd.concat([df_RNA_res,df_BER_res],axis=1)
df_RNA_res=pd.merge(df_RNA_res,res_pair1,on=['chrom','start'],how='left')
df_RNA_res.rename(columns={'score':'RNA_pair_st'},inplace=True)
df_RNA_res=pd.merge(df_RNA_res,res_pair2,on=['chrom','end'],how='left')
df_RNA_res.rename(columns={'score':'RNA_pair_ed'},inplace=True)
df_RNA_res=pd.merge(df_RNA_res,res_FE1,on=['chrom','start'],how='left')
df_RNA_res.rename(columns={'energ':'RNA_FE_st'},inplace=True)
df_RNA_res=pd.merge(df_RNA_res,res_FE2,on=['chrom','end'],how='left')
df_RNA_res.rename(columns={'energ':'RNA_FE_ed'},inplace=True)
fn_output=sys.argv[2]
df_RNA_res.to_csv(fn_output,sep='\t',index=False)

#initiate
os.system('rm {}* -rf'.format(dir_out))
