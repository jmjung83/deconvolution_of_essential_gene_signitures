'''
Created on 2018. 7. 23.

@author: jmjung
'''
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.kde import gaussian_kde
from numpy import linspace

if __name__ == '__main__':
    ## EG score

    eg_files = glob.glob('result/spearman/all_*_24H.txt')
    
    eg_df_list=[]
    eg_pval_df_list=[]
    for eg_file in eg_files:
        if 'merged' in eg_file:
            continue
        cell = eg_file.split('\\')[1].split('_')[1]
        eg_df = pd.read_table(eg_file, sep='\t', index_col = 0)
        eg_df['Essentiality Score (ES)'] = eg_df['EG_score'].round(3)
        eg_df.drop('EG_score',inplace=True, axis=1)
        eg_df = eg_df[['#_of_exp','#_of_exp(fc<0)','Essentiality Score (ES)','empirical_pVal']]
        eg_df.columns=eg_df.columns+" ({})".format(cell)
        
        eg_df_list.append(eg_df)
        
        eg_pval_df = eg_df['empirical_pVal ({})'.format(cell)].to_frame()
        eg_pval_df_list.append(eg_pval_df)
        
    merged_df=pd.concat(eg_df_list,axis=1)
    merged_pval_df=pd.concat(eg_pval_df_list,axis=1)
    merged_pval_df['Ave_pVal'] = merged_pval_df.mean(axis=1)
    merged_pval_df.sort_values(['Ave_pVal'], inplace=True)
    pd.DataFrame.to_csv(merged_df,path_or_buf='result/spearman/all_merged_24H.txt',sep='\t')
    pd.DataFrame.to_csv(merged_pval_df,path_or_buf='result/spearman/all_merged_24H_pValOnly.txt',sep='\t')

    pval_df = pd.read_table("result/spearman/all_merged_24H_pValOnly.txt", sep='\t', index_col = 0)
    for col in list(pval_df):
        print(col.replace("empirical_pVal",""), sum(pval_df[col]<=0.00001), sum(pval_df[col]<=0.00001)/12716.0)