'''
CorrieCaller v1.0.0
Juan C. Rueda-Silva, 2026
PARAMETERS
output_dir= name of the directory to save the output
test= name of the test group (should appear in the DiffBind output after one of the Conc_ columns)
control= name of the control group (should appear in the DiffBind output after one of the Conc_ columns)
diffBind_path= path to the DiffBind output file
sWindow_size= Size of the sliding window in bp (default: 10000)
sWindow_shift= Size of the sift bewteen sliding windows (default: 100)
modification_baseLine= Minum average WT (control) H3K9me3 enrichment per bin (default: 2)
bins_toCompare= the number of upstream/downstream bins to compare (default: 250)
cpus=Number of CPUs to use (default: 2)
'''

###
#INSERT your parameters bellow
###
output_dir='' #Directory to save the output (string)
test='' #Identifier of the "test" group (string)
control='' #Identifier of the "control" group (string)
diffBind_path='' #Path to the DiffBind output (string)
sWindow_size=10000 #Sliding window size in bp (integer)
sWindow_shift=100 #Difference betwween the start of subsequent sliding windows in bp (integer)
modification_baseLine=2 #Threshold for a bin to be considered H3K9me3-rich (float)
bins_toCompare=250 #Number of bins upstream and downstream to compare to call potential local maxima and minima (integer)
cpus=2 #Number of CPUs to use (integer)
###

import sys
import os
import pandas as pd
import multiprocessing as mp
from tqdm import tqdm
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def main():
    os.system('mkdir .//'+output_dir)
    path='.//'+output_dir+'//'
    global df_diffBind
    df_diffBind=pd.read_csv(diffBind_path,sep='\t',index_col=None)
    df_diffBind['chr']=df_diffBind['seqnames']
    global chrs
    chrs=list(set(df_diffBind['chr']))
    chrs_sizes=[max(list(df_diffBind.loc[df_diffBind['chr']==i,'end'])) for i in chrs]
    global chrs_sizes_map
    chrs_sizes_map=dict(zip(chrs,chrs_sizes))
    if cpus>4:
        processes=cpus-2
    else:
        processes=1
    os.system('echo "\nGenerating sliding window"')
    global sWindow_bins
    pool=mp.Pool(processes=processes)
    sWindow_bins=get_slidingWindow(chrs,chrs_sizes_map,sWindow_size,sWindow_shift,pool)
    pool.close()
    global sWindow_val
    os.system('echo "\nCalculating sliding window bins"')
    pool=mp.Pool(processes=processes)
    sWindow_val=get_sWindow_values(pool)
    pool.close()
    sWindow_val.to_csv(path+'SlidingWindow.tsv',sep='\t',index=False)
    global local_max_min_pred
    os.system('echo "\nEvaluating sliding window bins"')
    pool=mp.Pool(processes=processes)
    local_max_min_pred=predict_local_maxima_minima(pool)
    pool.close()
    global point_max_min
    os.system('echo "\nEstimmating local maxima and minima"')
    pool=mp.Pool(processes=processes)
    point_max_min=get_ponit_maxima_minima(pool)
    pool.close()
    global slopes
    os.system('echo "\nEstimmating slopes"')
    pool=mp.Pool(processes=processes)
    slopes=get_slopes(pool)
    pool.close()
    os.system('echo "\nCalling corries"')
    pool=mp.Pool(processes=processes)
    corries=call_corries(pool)
    pool.close()
    corries.to_csv(path+'corries.tsv',sep='\t',index=False)
    global reg_max_min_extra
    os.system('echo "\nReporting regions predicted to be local maxima/minima"')
    pool=mp.Pool(processes=processes)
    reg_max_min_extra=get_reg_max_min(pool)
    pool.close()
    reg_max_min=reg_max_min_extra.loc[(reg_max_min_extra['type']=='maximum')|(reg_max_min_extra['type']=='minimum')].copy()
    reg_max_min.to_csv(path+'predicted_local_maxima_minima.tsv',sep='\t',index=False)
    os.system('echo "\nReporting enriched regions"')
    enriched_regs=get_enriched_regs(reg_max_min_extra)
    enriched_regs.to_csv(path+'enriched_regs.tsv',sep='\t',index=False)
    os.system('echo "\nDone!"')

def get_slidingWindow(chrs,chrs_sizes_map,sWindow_size,sWindow_shift,pool):
    res=list(tqdm(pool.imap_unordered(get_sWindow_perChr,chrs),total=len(chrs)))
    df_res=pd.DataFrame(columns=['chr','start','end'])
    df_res=merge_res_df(df_res,res,list(df_res.columns))
    df_res=df_res.sort_values(by=['chr','start'])
    df_res=df_res.reset_index(drop=True)
    return df_res

def merge_res_df(df_res,res,columns):
    for r in res:
        sub_df=pd.DataFrame()
        for i in range(len(r)):
            sub_df[columns[i]]=r[i]
        df_res=pd.concat(\
                         [\
                          df_res,\
                          sub_df\
                         ],ignore_index=True\
                        )
    return df_res

def get_sWindow_perChr(c):
    starts=list(range(0,chrs_sizes_map[c],sWindow_shift))
    ends=[i+sWindow_size-1 for i in starts]
    return [[c]*len(starts),starts,ends]


def get_sWindow_values(pool):
    res=list(tqdm(pool.imap_unordered(get_sWindow_binVal,list(sWindow_bins.index.values)),total=len(sWindow_bins.index)))
    df_res=pd.DataFrame(res,columns=['chr','start','end','avg_l2fc','std_l2fc','avg_control','std_control'])
    df_res=df_res.sort_values(by=['chr','start'])
    df_res=df_res.reset_index(drop=True)
    return df_res

def get_sWindow_binVal(i):
    chrom=sWindow_bins.loc[i,'chr']
    start=sWindow_bins.loc[i,'start']
    end=sWindow_bins.loc[i,'end']
    sub_df=get_intersect_df(df_diffBind,chrom,start,end)
    sub_df.loc[sub_df['start']<start,'start']=start
    sub_df.loc[sub_df['end']>end,'end']=end
    sub_df['width']=sub_df['end']-sub_df['start']
    sub_df['weights']=sub_df['width']/(end-start)
    avg_fold=get_avg(sub_df,'Fold')
    std_fold=get_std(sub_df,'Fold',avg_fold)
    avg_c=get_avg(sub_df,'Conc_'+control)
    std_c=get_std(sub_df,'Conc_'+control,avg_c)
    return [chrom,start,end,avg_fold,std_fold,avg_c,std_c]

def get_avg(df,col):
    return np.sum(df['weights']*df[col])

def get_std(df,col,avg):
    return np.sqrt(np.sum(df['weights']*np.power((df[col]-avg),2))/(len(df[col]-1)*np.sum(df['weights'])/len(df[col])))

def get_intersect_df(df,chrom,start,end):
    df=df.loc[df['chr']==chrom].copy()
    sub_df=df.loc[((df['start']>=start) & (df['end']<=end)) | \
                  ((df['start']<=start) & (df['end']>=start)) | \
                  ((df['start']<=end) & (df['end']>=end))|\
                  ((df['start']<=start) & (df['end']>=start))|\
                  ((df['start']>=start) & (df['start']<=end))|\
                  ((df['end']>=start) & (df['end']<=end))\
                 ].copy()
    return sub_df

def predict_local_maxima_minima(pool):
    res=list(tqdm(pool.imap_unordered(get_bin_evaluation,chrs),total=len(chrs)))
    df_res=pd.DataFrame(columns=['chr','start','end','evaluation'])
    df_res=merge_res_df(df_res,res,list(df_res.columns))
    df_res=df_res.sort_values(by=['chr','start'])
    return df_res

def get_bin_evaluation(c):
    df=sWindow_val.loc[sWindow_val['chr']==c].copy()
    wt=list(df['avg_control'])
    df['evaluation']=get_max(list(df['avg_l2fc']),wt)
    return [list(df['chr']),list(df['start']),list(df['end']),list(df['evaluation'])]

def get_max(l,wt):
    res=['outside']*bins_toCompare
    for i in range(bins_toCompare,len(l)-bins_toCompare):
        res.append(get_type(l,i,wt))
    return res+['outside']*bins_toCompare

def get_type(l,i,wt):
    b=l[i]
    if wt[i]>modification_baseLine:
        up=np.mean(l[i-bins_toCompare:i])
        down=np.mean(l[i+1:i+bins_toCompare])
        if b>up and b>down:
            return 'maximum'
        elif b<down and b<up:
            return 'minimum'
        else:
            return 'none'
    else:
        return 'outside'
def get_ponit_maxima_minima(pool):
    res=list(tqdm(pool.imap_unordered(get_point_eval,chrs),total=len(chrs)))
    df_res=pd.DataFrame(columns=['chr','start','pos','end','type'])
    df_res=merge_res_df(df_res,res,list(df_res.columns))
    df_res=df_res.sort_values(by=['chr','start'])
    df_res=df_res.reset_index(drop=True)
    return df_res

def get_point_eval(chrom):
    df=local_max_min_pred.loc[local_max_min_pred['chr']==chrom].copy()
    df=df.reset_index()
    df['midpoint']=df['start']+(df['end']-df['start'])/2
    df_res=pd.DataFrame(get_vals_point(df,chrom),columns=['chr','start','pos','end','type'])
    return [list(df_res['chr']),list(df_res['start']),list(df_res['pos']),list(df_res['end']),list(df_res['type'])]

def get_vals_point(df,chrom):
    c='evaluation'
    out=True
    vals=[]
    res=[]
    row=[]
    for i in list(df.index.values):
        vals.append(df.loc[i,'midpoint'])
        if i!=len(df.index)-1:
            if out:
                out=False
                if df.loc[i,c]!=df.loc[i+1,c]:
                    row.append(chrom)
                    row.append(int(min(vals)))
                    row.append(int(np.mean(vals)))
                    row.append(int(max(vals)))
                    row.append(df.loc[i,c])
                    res.append(row)
                    out=True
                    vals=[]
                    row=[]
            else:
                if df.loc[i,c]!=df.loc[i+1,c]:
                    row.append(chrom)
                    row.append(int(min(vals)))
                    row.append(int(np.mean(vals)))
                    row.append(int(max(vals)))
                    row.append(df.loc[i,c])
                    res.append(row)
                    out=True
                    vals=[]
                    row=[]
        else:
            row.append(chrom)
            row.append(int(min(vals)))
            row.append(int(np.mean(vals)))
            row.append(int(max(vals)))
            row.append(df.loc[i,c])
            res.append(row)
            out=True
            vals=[]
            row=[]
    return res

def get_slopes(pool):
    res=list(tqdm(pool.imap_unordered(slopes_calc,chrs),total=len(chrs)))
    df_res=pd.DataFrame(columns=['chr','start_extended','start','start_strict','end_strict','end','end_extended','strand','type'])
    df_res=merge_res_df(df_res,res,list(df_res.columns))
    df_res=df_res.sort_values(by=['chr','start'])
    df_res=df_res.reset_index(drop=True)
    return df_res

def slopes_calc(chrom):
    df=point_max_min.loc[point_max_min['chr']==chrom].copy()
    df_res=pd.DataFrame(find_slopes(df),columns=['chr','start_extended','start','start_strict','end_strict','end','end_extended','strand','type'])
    return [list(df_res[i]) for i in list(df_res.columns)]

def find_slopes(df):
    window=sWindow_size
    inReg=False
    res=[]
    row=[]
    current=''
    for i in list(df.index.values):
        if inReg:
            if i==len(df.index) or df.loc[i,'type']=='outside':
                if df.loc[i,'type']!='none' and df.loc[i,'type']!='outside':
                    type='slope'
                    if current=='minimum':
                        strand='+'
                    else:
                        strand='-'
                else:
                    strand='*'
                    if current=='minimum':
                        type='minimum'
                    else:
                        type='maximum'
                row.append(df.loc[i,'start'])
                row.append(df.loc[i,'pos'])
                row.append(int(df.loc[i,'end']+(window/2)))
                row.append(strand)
                row.append(type)
                res.append(row)
                row=[]
                inReg=False
                current=''
            elif df.loc[i,'type']!='none' and df.loc[i,'type']!=current:
                type='slope'
                if current=='minimum' and df.loc[i,'type']=='maximum':
                    strand='+'
                else:
                    strand='-'
                row.append(df.loc[i,'start'])
                row.append(df.loc[i,'pos'])
                row.append(int(df.loc[i,'end']+(window/2)))
                row.append(strand)
                row.append(type)
                res.append(row)
                row=[]
                inReg=True
                current=df.loc[i,'type']
                row.append(df.loc[i,'chr'])
                row.append(int(df.loc[i,'start']-(window/2)))
                row.append(df.loc[i,'pos'])
                row.append(df.loc[i,'end'])
            elif df.loc[i,'type']!='none' and df.loc[i,'type']==current:
                strand='*'
                if current=='minimum':
                    type='minimum'
                else:
                    type='maximum'
                row.append(df.loc[i,'start'])
                row.append(df.loc[i,'pos'])
                row.append(int(df.loc[i,'end']+(window/2)))
                row.append(strand)
                row.append(type)
                res.append(row)
                row=[]
                inReg=True
                current=df.loc[i,'type']
                row.append(df.loc[i,'chr'])
                row.append(int(df.loc[i,'start']-(window/2)))
                row.append(df.loc[i,'pos'])
                row.append(df.loc[i,'end'])

        else:
            if df.loc[i,'type']!='none' and df.loc[i,'type']!='outside':
                inReg=True
                current=df.loc[i,'type']
                row.append(df.loc[i,'chr'])
                row.append(int(df.loc[i,'start']-(window/2)))
                row.append(df.loc[i,'pos'])
                row.append(df.loc[i,'end'])
    return res

def call_corries(pool):
    res=list(tqdm(pool.imap_unordered(get_corrie,list(slopes.index.values)),total=len(slopes.index)))
    res=[i for i in res if i!=False]
    df_res=pd.DataFrame(res,columns=['chr','start_extended','start','start_strict','end_strict','end','end_extended','strand'])
    df_res['width']=df_res['end']-df_res['start']
    df_res=df_res.sort_values(by=['chr','start'])
    df_res=df_res.reset_index(drop=True)
    return df_res

def get_corrie(i):
    indexes=list(slopes.index.values)
    df=slopes
    if i!=indexes[0] and i!=indexes[-1]:
        chrom=df.loc[i,'chr']
        if chrom==df.loc[i+1,'chr'] and chrom==df.loc[i-1,'chr']:
            type=df.loc[i,'type']
            if type=='slope':
                return [chrom,df.loc[i,'start_extended'],df.loc[i,'start'],df.loc[i,'start_strict'],df.loc[i,'end_strict'],\
                        df.loc[i,'end'],df.loc[i,'end_extended'],df.loc[i,'strand']]
            else:
                return False
        else:
            return False
    else:
        return False

def get_reg_max_min(pool):
    res=list(tqdm(pool.imap_unordered(cacl_reg_min_max,chrs),total=len(chrs)))
    df_res=pd.DataFrame(columns=['chr','start','end','type'])
    df_res=merge_res_df(df_res,res,list(df_res.columns))
    df_res=df_res.sort_values(by=['chr','start'])
    df_res=df_res.reset_index(drop=True)
    return df_res

def cacl_reg_min_max(chrom):
    df=local_max_min_pred.loc[local_max_min_pred['chr']==chrom].copy()
    df['midpoint']=df['start']+(df['end']-df['start'])/2
    df_res=pd.DataFrame(get_vals_regs(df,chrom),columns=['chr','start','end','type'])
    return [list(df_res[i]) for i in list(df_res.columns)]

def get_vals_regs(df,chrom):
    c='evaluation'
    df=df.reset_index()
    out=True
    vals=[]
    res=[]
    row=[]
    for i in list(df.index.values):
        if i!=len(df.index)-1:
            if out:
                row.append(chrom)
                row.append(int(df.loc[i,'midpoint']))
                out=False
                if df.loc[i,c]!=df.loc[i+1,c]:
                    row.append(int(df.loc[i,'midpoint']))
                    row.append(df.loc[i,c])
                    res.append(row)
                    out=True
                    vals=[]
                    row=[]
            else:
                if df.loc[i,c]!=df.loc[i+1,c]:
                    row.append(int(df.loc[i,'midpoint']))
                    row.append(df.loc[i,c])
                    res.append(row)
                    out=True
                    vals=[]
                    row=[]
        else:
            row.append(int(df.loc[i,'midpoint']))
            row.append(df.loc[i,c])
            res.append(row)
            out=True
            vals=[]
            row=[]
    return res

def get_enriched_regs(reg_max_min_extra):
    df=reg_max_min_extra.loc[reg_max_min_extra['type']!='outside'].copy()
    inReg=False
    chrom=''
    start=0
    end=0
    res=[]
    indexes=list(df.index.values)

    for n in range(len(indexes)-1):
        i=indexes[n]
        i_n=indexes[n+1]
        if not(inReg):
            chrom=df.loc[i,'chr']
            start=df.loc[i,'start']
            inReg=True
        if inReg:
            if chrom!=df.loc[i_n,'chr'] or df.loc[i,'end']!=df.loc[i_n,'start']-100:
                res.append([chrom,start,df.loc[i,'end']])
                inReg=False
    if inReg:
        res.append([chrom,start,df.loc[indexes[-1],'end']])
    else:
        res.append([df.loc[indexes[-1],'chr'],df.loc[indexes[-1],'start'],df.loc[indexes[-1],'end']])

    df_res=pd.DataFrame(res,columns=['chr','start','end'])
    df_res['width']=df_res['end']-df_res['start']
    df_res=df_res.sort_values(by=['chr','start'])
    df_res=df_res.reset_index(drop=True)
    return df_res

if __name__=='__main__':
    main()
