'''
Created on Feb 24, 2017

@author: jmjung
'''
import pandas as pd
import time
from commonFunc import getCols
import numpy as np
import glob
def measure(init_time):    
    after_time=time.time()
    dif_time=after_time-init_time                                     
    hour=int(dif_time/3660)
    mins=int((dif_time-hour*3660)/60)
    sec=dif_time-hour*3660-mins*60                                 
    print 'Processing Time:'+str(hour) +' hour\t'+str(mins) +' min\t'+str(sec) +' sec'

def quantileNormalize(df_input):
#     #another algorithm
#     rank_mean = subLns_df.stack().groupby(subLns_df.rank(method='first').stack().astype(int)).mean()
#     qtSubLns_df = subLns_df.rank(method='min').stack().astype(int).map(rank_mean).unstack()

    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

def getDirectProbe():
    probe_file = open('../rawData/probe_directed_GPL20573.txt')

    directProbe_list = []
    for line in probe_file.xreadlines():
        probe = line.replace('\r','').replace('\n','')
        directProbe_list.append(probe)
        
    probe_file.close()
    return directProbe_list

if __name__ == '__main__':
    all_files = glob.glob('preProcessResult/LINCS/*.txt')
    for in_file in all_files:
        if in_file.split('_')[1] != '24':
            continue
        for mode in ['all', 'lmk']:
            print in_file[23:-4], mode
            init_time = time.time()
            
            # read expression data in LINCS per cell-hour
            lns_df = pd.read_table(in_file, sep='\t', index_col = 0, dtype = 'unicode')
            
            # get distinct plate name
            pt_list = [pt for pt in list(set(lns_df.ix['det_plate'])) if pt != '-666']
            
            compDoseFC_dic = {}
            compDoseFCcnt_dic = {}
            for pt in pt_list:
                # get only expression data
                lnsDat_df = lns_df.ix[12:,lns_df.ix['det_plate']==pt]
                # convert data type to numeric
                lnsDat_df = lnsDat_df.convert_objects(convert_numeric=True)
                # quantileNormalization per plate
                lnsDat_df = quantileNormalize(lnsDat_df)
                
                # attach gene info
                lnsG_df = lns_df.ix[12:,:2]
                lnsGdat_df = pd.concat([lnsG_df, lnsDat_df], axis=1)
                #print lnsDat_df.ix[:20,:20].to_string()
                
                # get expression average per gene for two mode, all vs landmark gene
                if mode == 'all':
                    lnsAveDat_df = lnsGdat_df.groupby('pr_gene_symbol').mean()
                    lnsAveDat_df.drop('-666', inplace=True)
                elif mode == 'lmk':
                    directProbe_list = getDirectProbe()
                    lnsAveDat_df = lnsGdat_df.ix[directProbe_list].groupby('pr_gene_symbol').mean()
                
                # attach experiment info
                lnsInfo_df = lns_df.ix[:12,lns_df.ix['det_plate']==pt]
                lnsAve_df = pd.concat([lnsInfo_df,lnsAveDat_df])
                
                #print lnsAve_df.ix[:20,:20].to_string()
                
                # get gene expression average of all DMSO
                lnsAveDMSO_df = lnsAve_df.ix[:,(lnsAve_df.ix['sm_lincs_id']=='DMSO')|(lnsAve_df.ix['sm_lincs_id']=='dmso')]
                lnsAveDMSOmean_ds = lnsAveDMSO_df.ix[12:].mean(axis=1)
                
                lnsAveTreat_df= lnsAve_df.ix[:,(lnsAve_df.ix['sm_lincs_id']!='DMSO')&(lnsAve_df.ix['sm_lincs_id']!='dmso')]
                
                
                for exp in list(lnsAveTreat_df):
                    
                    comp = lnsAveTreat_df.ix['sm_lincs_id',exp]
                    dose = lnsAveTreat_df.ix['sm_dose',exp]
                    
                    # get fold change per each experiment, values are already log scale
                    FC_ds = lnsAveTreat_df.ix[12:,exp] - lnsAveDMSOmean_ds
                    
                    if FC_ds.isnull().values.any():
                        print '############## NAN exist', in_file, pt, exp
                    
                    # get mean fold change per comp-dose pair
                    if not compDoseFC_dic.has_key((comp, dose)):
                        compDoseFC_dic[(comp, dose)] = FC_ds
                        compDoseFCcnt_dic[(comp, dose)] = 1
                    else:
                        compDoseFC_dic[(comp, dose)] += FC_ds
                        compDoseFCcnt_dic[(comp, dose)] += 1
            
            compDoseFCave_dic = {}
            for key in compDoseFC_dic.keys():
                compDoseFCave_dic['|'.join(key)] = compDoseFC_dic[key]/float(compDoseFCcnt_dic[key])
                
            compDose_df = pd.DataFrame(compDoseFCave_dic)
            ## write to tsv file
            pd.DataFrame.to_csv(compDose_df, path_or_buf='result/foldChange/%s_%sH.txt'%(mode, in_file[23:-10]),sep='\t')
     
            measure(init_time)
        