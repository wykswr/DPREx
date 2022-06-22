import pandas as pd

W=100 #window size

#chromosome size
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

#load bed_file_ori
#for check code
#fn_bed1='./preprocess/gnomAD-repeats_pos.bed'

#get_newposition
def new_pos(df_bed_ori,nr): #nr--number of row
    chrr,st,ed=df_bed_ori.loc[nr]
    edmax=HG19_CHROMSIZE[chrr]
    pos1=max([0,st-W])
    pos2=min([edmax,ed+W])
    return([(pos1,st),(ed,pos2)])

def new_pos1(df_bed_ori):
    l_pos1=[]
    l_pos2=[]
    for nr in df_bed_ori.index:
        [(pos1,_),(_,pos2)]=new_pos(df_bed_ori,nr)
        l_pos1.append(pos1)
        l_pos2.append(pos2)
    df1=pd.DataFrame({'chr':df_bed_ori.chrom,'st':l_pos1,'ed':df_bed_ori.start})
    df2=pd.DataFrame({'chr':df_bed_ori.chrom,'st':df_bed_ori.end,'ed':l_pos2})
    return(df1,df2)

