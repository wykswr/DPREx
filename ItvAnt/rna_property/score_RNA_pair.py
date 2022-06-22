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
    df_score=pd.DataFrame({'chrom':l_chr,'pos':l_pos,'score':l_score})
    df_score['chrom']=df_score.chrom.str.lower()
    df_score[['pos1','pos2']]=df_score.pos.str.split('-',expand=True)
    ncol=['pos2','pos1'][n]
    df_score=df_score[['chrom',ncol,'score']]
    df_score[ncol]=df_score[ncol].astype('int64')
    df_score.rename(columns={ncol:['start','end'][n]},inplace=True)
    return(df_score)


