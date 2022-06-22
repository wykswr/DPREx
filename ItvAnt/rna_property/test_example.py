#!/home/fancong/anaconda3/bin/python3.8
import pandas as pd
import main

#read input
def read_input(fn_in:str):
    pd_loc=pd.read_csv(fn_in,header=None,names=['chrom','start','end'],sep='\t')
    pd_loc.iloc[:,[1,2]]=pd_loc.iloc[:,[1,2]].astype('int64')
    return(pd_loc)

pd_in=read_input('./example.tsv')
pd_res=main.get_rna_property(pd_in,'/bigdat1/user/fancong/RNA_charac_cal')

pd_res.to_csv('fanc_test_res.tsv',sep='\t',index=False)
