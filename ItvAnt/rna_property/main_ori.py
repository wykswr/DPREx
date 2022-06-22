#!/home/fancong/anaconda3/bin/python3.8
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
from argparse import ArgumentParser

#calculation paths 
with open(Path(__file__).parent / 'fan_path11.json') as handle:
    config = json.load(handle)
bash=config["bash"]
bedtools_run=config["bedtools"]
conda_python3=config['conda_python3']
#accept parameters
def get_args():
    parser = ArgumentParser(description="I've been laboring on this damn project for more than one year")
    parser.add_argument('--input_file', help='in hg19 reference gnome, chr, start, end of fragment, separated by tab')
    parser.add_argument('--output_file', help='chr, start, end, and values of index, separated by tab')
    parser.add_argument('--dir_rloop_annot', help='path of the hg19.rloop.filtered.mergedRegion.bed')
    parser.add_argument('--dir_BER_annot', help='path of the four files of ACNEIL_SSINP_2_hg19_FE,GSM2137770_POLB_hg19,ogg1_hg19 and xrcc1_hg19')
    parser.add_argument('--hg19_in', help='the position of the hg19 reference genome in fasta format, NOTE, the corresponding .fai file should be in the same path')
    parser.add_argument('--dir_GF',help='path you install GF_GRASS, the path contains filles like train.py, test.py and utils.py, etc.')
    parser.add_argument('--run_Rfold', help='the executive script of Rloop')
    parser.add_argument('--dir_gen_files', help='the path to contain the generated files during calculation')
    return parser.parse_args()

#run shell script
def run_bash(cmd):
    p = Popen([bash, '-c', cmd], stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    out, err = out.decode('utf8'), err.decode('utf8')
    rc = p.returncode
    if rc != 0:
        print(err)
    return rc, out, err

#read input
def read_input(fn_in:str):
    pd_loc=pd.read_csv(fn_in,header=None,names=['chrom','start','end'],sep='\t')
    pd_loc.iloc[:,[1,2]]=pd_loc.iloc[:,[1,2]].astype('int64')
    return(pd_loc)

def annot_Rloop_BER(pd_in,dir_rloop_annot,dir_BER_annot):
    #annot_Rloop
    fn_rloop_annot= Path(dir_rloop_annot) / "hg19.rloop.filtered.mergedRegion.bed"
    l_rloop=pd_in.index.map(lambda x:annot_Rloop.get_olp2(x,pd_in,fn_rloop_annot))
    #calculate BER, DNA_repair factors
    pd_rloop_BER=annot_BER.cal_sum1(pd_in,dir_BER_annot)
    #combine result
    pd_rloop_BER['Rloop']=l_rloop
    return(pd_rloop_BER)
    
def input2seq(pd_in,dir_gen_files,hg19_in): 
    #transfer to sequence
    df1,df2=get_bed.new_pos1(pd_in) #entend 100bp
    (Path(dir_gen_files) / "new_bed").mkdir(exist_ok=True)
    dir_bed=Path(dir_gen_files) / "new_bed"
    df1.to_csv(Path(dir_bed) / 'df1.bed',index=False,header=False,sep='\t')
    df2.to_csv(Path(dir_bed) / 'df2.bed',index=False,header=False,sep='\t')
    (Path(dir_gen_files) / 'get_seq').mkdir(exist_ok=True)
    dir_seq=Path(dir_gen_files) / 'get_seq'
    run_bash('{} getfasta -fi {} -bed {} -fo {}'.format(bedtools_run,hg19_in,(Path(dir_bed) / 'df1.bed'),Path(dir_seq) / 'df1.fa'))
    run_bash('{} getfasta -fi {} -bed {} -fo {}'.format(bedtools_run,hg19_in,(Path(dir_bed) / 'df2.bed'),Path(dir_seq) / 'df2.fa'))
    return()

def cal_RNA_pair(dir_gen_files,dir_GF,conda_python3):
    (Path(dir_gen_files) / 'RNA_pair').mkdir(exist_ok=True) 
    dir_RNA_pair=Path(dir_gen_files) / 'RNA_pair'
    dir_seq=Path(dir_gen_files) / 'get_seq'
    run_bash('{} {} {} {} {}'.format('./run_GRASP.sh',dir_GF,dir_seq,dir_RNA_pair,conda_python3))
    fn_df1=Path(dir_RNA_pair) / 'df1_upper.out'
    fn_df2=Path(dir_RNA_pair) / 'df2_upper.out'
    res_pair1=score_RNA_pair.read_score(fn_df1,0) #start
    res_pair2=score_RNA_pair.read_score(fn_df2,1) #end
    return(res_pair1,res_pair2)

def cal_RNA_FE(dir_gen_files,run_RNA_fold):
    dir_seq=Path(dir_gen_files) / 'get_seq'
    (Path(dir_gen_files) / 'RNA_free_energy').mkdir(exist_ok=True)
    dir_RNA_FE=Path(dir_gen_files) / 'RNA_free_energy'
    fn_seq1=Path(dir_seq) / 'df1.fa'
    fn_seq2=Path(dir_seq) / 'df2.fa'
    run_bash('{} {} {} {}'.format('./RNA_FE.sh',fn_seq1,dir_RNA_FE,run_RNA_fold))
    run_bash('{} {} {} {}'.format('./RNA_FE.sh',fn_seq2,dir_RNA_FE,run_RNA_fold))
    res_FE1=pd.read_csv(Path(dir_RNA_FE / 'df1.energ'),sep='\t|-|_',header=None,names=['chrom','start','energ'],usecols=[0,2,4],engine='python')
    res_FE2=pd.read_csv(Path(dir_RNA_FE / 'df2.energ'),sep='\t|-|_',header=None,names=['chrom','end','energ'],usecols=[0,1,4],engine='python')
    return(res_FE1,res_FE2)

def annot_RNA():
    pd_in=read_input(fn_in)
    pd_rloop_BER=annot_Rloop_BER(pd_in,dir_rloop_annot,dir_BER_annot)
    input2seq(pd_in,dir_gen_files,hg19_in)
    res_pair1,res_pair2=cal_RNA_pair(dir_gen_files,dir_GF,conda_python3)
    res_FE1,res_FE2=cal_RNA_FE(dir_gen_files,run_Rfold)
    #combine result
    df_RNA_res=pd_in.loc[:,:]
    df_RNA_res=pd.concat([df_RNA_res,pd_rloop_BER],axis=1)
    df_RNA_res=pd.merge(df_RNA_res,res_pair1,on=['chrom','start'],how='left')
    df_RNA_res.rename(columns={'score':'RNA_pair_st'},inplace=True)
    df_RNA_res=pd.merge(df_RNA_res,res_pair2,on=['chrom','end'],how='left')
    df_RNA_res.rename(columns={'score':'RNA_pair_ed'},inplace=True)
    df_RNA_res=pd.merge(df_RNA_res,res_FE1,on=['chrom','start'],how='left')
    df_RNA_res.rename(columns={'energ':'RNA_FE_st'},inplace=True)
    df_RNA_res=pd.merge(df_RNA_res,res_FE2,on=['chrom','end'],how='left')
    df_RNA_res.rename(columns={'energ':'RNA_FE_ed'},inplace=True)
    df_RNA_res.to_csv(fn_out,sep='\t',index=False)
    return()

def deleted_gen_file():
    os.system('rm {}* -rf'.format(dir_gen_files))
    return()

if __name__ == '__main__':
    args = get_args()
    fn_in=args.input_file
    fn_out=args.output_file
    dir_rloop_annot=args.dir_rloop_annot
    dir_BER_annot=args.dir_BER_annot
    hg19_in=args.hg19_in
    dir_GF=args.dir_GF
    run_Rfold=args.run_Rfold
    dir_gen_files=args.dir_gen_files
    annot_RNA()
    #deleted_gen_file()

