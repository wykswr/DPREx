import pandas as pd
import pathlib
from pathlib import Path

#calculated overlap percentage
def read_rloop_annot(fn_rloop_annot):
    df_rloop=pd.read_csv(fn_rloop_annot,sep='\t',usecols=[0,1,2],header=None,names=['chr','st','ed'])
    return(df_rloop)

def get_olp1(st1,ed1,st2,ed2): #1--gfr,2--repeat
    #the range had been assigned in other steps!!!
    lrep=ed2-st2 #length of repeat,note!length from 121-122 is 2bp!
    lolp=0 #default length of overlap
    l=[st1,st2,ed1,ed2]
    l.sort()
    lolp=l[2]-l[1] #true length of overlap
    p=lolp/lrep
    return(p)

def get_olp2(ind,df_bed,fn_rloop_annot): #ind--index
    df_rloop=read_rloop_annot(fn_rloop_annot)
    #set the range of st1, ed1, st2 and ed2 here!
    l_p=[0]
    nchr,st2,ed2 = df_bed.iloc[ind,[0,1,2]]
    df=df_rloop[df_rloop.chr==nchr]
    df=df[~((df.ed <= st2)|(df.st >= ed2))]
    if df.shape[0]>0:
        l_p=df.index.map(lambda x: get_olp1(df.loc[x,'st'],df.loc[x,'ed'],st2,ed2))
    return(max(l_p))
