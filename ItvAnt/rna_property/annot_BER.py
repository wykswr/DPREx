import pyBigWig
import numpy as np
import pandas as pd
import pathlib
from pathlib import Path

HG19_CHROMSIZE = {"chr1": 249250621, "chr2": 243199373,
        "chr3": 198022430, "chr4": 191154276,
        "chr5": 180915260, "chr6": 171115067,
        "chr7": 159138663, "chr8": 146364022,
        "chr9": 141213431, "chr10": 135534747,
        "chr11": 135006516, "chr12": 133851895,
        "chr13": 115169878, "chr14": 107349540,
        "chr15": 102531392, "chr16": 90354753,
        "chr17": 81195210, "chr18": 78077248,
        "chr19": 59128983, "chr20": 63025520,
        "chr21": 48129895, "chr22": 51304566,
        "chrX": 155270560, "chrY": 59373566,
        "chrM": 16569, "chrMT": 16569}

l_DNAr_bn=['ACNEIL_SSINP_2_hg19_FE','GSM2137770_POLB_hg19','ogg1_hg19','xrcc1_hg19']
#u=np.log10(3/10) #change!!! Strange!u should be positive!
u=np.log(19)
w=0.5

def cal_w(pos,ct):#pos--position of bigwig; ct--center. start or end
    d=abs(pos-ct)
    wd=2*np.exp(-u*d)/(1+np.exp(-u*d))
    return(wd)
    
def sum_signal(bw,indr,df,w): #w--window, 0.5kb, 1kb; indr--index of rows.
    nchr,st,ed=df.loc[indr,['chrom','start','end']]
    edmax=HG19_CHROMSIZE[nchr]
    #M=math.ceil(np.mean([st,ed]))
    st1=max(0,int(st-w*1000))
    st2=min(edmax,int(st+w*1000))
    ed1=max(0,ed-int(w*1000))
    ed2=min(edmax,int(ed+w*1000))
    l_value1=bw.values(nchr,st1,st2+1)
    l_w1=list(map(lambda x:cal_w(x,st),range(st1,st2+1)))
    l1=[x*y for x,y in zip(l_value1,l_w1) if ~(np.isnan(x))]
    l_value2=bw.values(nchr,ed1,ed2+1)
    l_w2=list(map(lambda x:cal_w(x,ed),range(ed1,ed2+1)))
    l2=[x*y for x,y in zip(l_value2,l_w2) if ~(np.isnan(x))]
    #print(w,len(l1),len(l2)) #length is almost the same as w, correct.
    sm=np.sum(l1)+np.sum(l2)
    return(sm)

def cal_sum1(df_bed,dir_BER_annot):
    ind_bw=0
    df_res=pd.DataFrame()
    while ind_bw<len(l_DNAr_bn):
        bn2=l_DNAr_bn[ind_bw]+'.bw'
        fn2=Path(dir_BER_annot) / bn2
        #print(str(fn2))
        bw=pyBigWig.open(str(fn2))
        l_sum=df_bed.index.map(lambda x:sum_signal(bw,x,df_bed,w))
        df_res['DNA_repair_'+str(ind_bw)]=l_sum
        bw.close()
        ind_bw+=1
    df_res=df_res.round(4)
    return(df_res)




