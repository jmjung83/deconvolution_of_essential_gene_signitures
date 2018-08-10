'''
Created on Sep 1, 2015

@author: jmjung
'''

import numpy as np
import math
import scipy.stats
import matplotlib.pyplot as plt
import pandas as pd
import glob

def getAvePoint(ge, gi):
    
    ITV=0.1
    giAve=[]
    geAve=[]
    geDisc = np.arange(int(min(ge))-1, int(max(ge))+1,ITV)
    #print(geDisc)
    for ii in geDisc:
        index_list = [ge.index(val) for val in ge if (val>=ii and val<(ii+0.2))]
        if len(index_list) == 0:
            continue
        else:
            geAve.append(ii+ITV/2.0)
            giAve.append(np.mean([gi[ix] for ix in index_list]))
    
    #print(geAve, giAve)
    return geAve, giAve

### Main Function ###
if __name__ == '__main__':
    sp_df = pd.read_table('result/spearman/all_A375_24H.txt', sep='\t', index_col = 0)
    pval_dic={}
    for gene in ['RPS19','TCOF1','NIP7','CARS2','RPL7']:
        pval_dic[gene] = [round(sp_df.loc[gene]['EG_score'],3),sp_df.loc[gene]['empirical_pVal'], int(sp_df.loc[gene]['#_of_exp(fc<0)'])]
    
    fc_df = pd.read_table('result/foldChange/all_A375_24H.txt', sep='\t', index_col = 0)
    cgi_df = pd.read_table('result/cellGrowthInh/A375_24H.txt', sep='\t', index_col = 0)
    
    plt.figure(figsize=(50,20))
    plotCnt=0
    for gene in ['RPS19','TCOF1','NIP7','CARS2','RPL7']:
        concat_df = pd.concat([fc_df, cgi_df])
        concat_df = concat_df.ix[:,concat_df.ix[gene]<0.0]
        ge_list = list(concat_df.ix[gene])
        gi_list = list(concat_df.ix['CGI'])
        print(gene, min(ge_list), max(gi_list))
        plotCnt += 1
        plt.subplot(2,5,plotCnt)
        plt.plot(ge_list, gi_list, 'o', mfc='none', markersize=8, color = 'dimgray')
        
        plt.title('A375_%s (negEXP: %s)\n(ES: %s, p-val:%s)'%(gene,pval_dic[gene][2],pval_dic[gene][0],pval_dic[gene][1]))
        
        plt.subplot(2,5,plotCnt + 5)
        geAve_list, giAve_list = getAvePoint(ge_list, gi_list)
        plt.plot(geAve_list, giAve_list, '-o', color = 'black')
        
    plt.show()    
    


        