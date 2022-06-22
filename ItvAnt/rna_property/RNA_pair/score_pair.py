#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 14:32:18 2021
process pair_prob,try mean value first
@author: fanc232
"""

import pandas as pd
creatVar=locals()

def read_score(fn,n:[0,1]): #"-32.4" in score contains "-" ; note! st-100 and ed+100
    l_chr=[]
    l_pos=[]
    l_score=[]
    with open (fn) as fm:
        for line in fm:
            if '>' in line:
                chrr,pos=line.strip().strip('>').split(':')
                l_chr.append(chrr)
                l_pos.append(pos)
            if ';' in line:
                l=line.strip().split(';')[0:100]
                l=[float(x) for x in l]
                score=sum(l)/100
                l_score.append(score)
    df_score=pd.DataFrame({'chr':l_chr,'pos':l_pos,'score':l_score})
    df_score.drop_duplicates(keep='first',inplace=True)
    df_score[['pos1','pos2']]=df_score.pos.str.split('-',expand=True)
    ncol=['pos2','pos1'][n]
    df_score=df_score[['chr',ncol,'score']]
    df_score[ncol]=df_score[ncol].astype('int64')
    df_score.rename(columns={ncol:['st','ed'][n]},inplace=True)
    fn_out=bn+'_pair_seq'+str(n)+'.score'
    df_score.to_csv(fn_out,index=False)
    return()

for bn in ['gnomAD','HGMD','NR','TRF']:
    fn1='/home/fanc232/tasks/ck_repeat/wangyuk/RNA_characters/test/RNA_pair/gf_grass/'+bn+'_for_seq1_upper.out'
    fn2='/home/fanc232/tasks/ck_repeat/wangyuk/RNA_characters/test/RNA_pair/gf_grass/'+bn+'_for_seq2_upper.out'
    read_score(fn1,0)
    read_score(fn2,1)

