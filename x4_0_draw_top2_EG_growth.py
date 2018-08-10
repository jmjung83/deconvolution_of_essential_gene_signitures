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



def plotContFC((g1, g1_list, g2, g2_list, gi_list, apr_file)):
    
    plotPos = 0
    for ge_list, ge in [(g1_list,g1), (g2_list,g2)]:
        plotPos += 1
        if plotPos == 3: plotPos=4
        plt.subplot(3,3,plotPos)
        plt.plot(ge_list, gi_list, 'o', mfc='none', markersize=5)
        plt.scatter(ge_list, [-20]*len(ge_list), s = 10, c=gi_list, cmap='Blues')
        
        plt.xlabel('fold change (%s)'%(ge), fontsize=12)
        plt.ylabel('growth inhibition', fontsize=12)
        
        plotPos += 1
        geAve_list, giAve_list = getAvePoint(ge_list, gi_list)
        plt.subplot(3,3,plotPos)
        plt.plot(geAve_list, giAve_list, '-o', mfc='none')
        plt.xlabel('fold change (%s)'%(ge), fontsize=12)
        plt.ylabel('growth Inh. (ave)', fontsize=12)

    plt.subplot(3,3,7)
    plt.scatter(g1_list, g2_list, s= 10, c=gi_list, cmap='Blues')
    plt.xlabel('fold change (%s)'%(g1), fontsize=12)
    plt.ylabel('fold change (%s)'%(g2), fontsize=12)
    
    g1Ave_list, g2Ave_list, giAve_list = get2dAvePoint(g1_list, g2_list, gi_list)
    plt.subplot(3,3,8)
    plt.scatter(g1Ave_list, g2Ave_list, s= 10, c=giAve_list, cmap='Blues')
    plt.xlabel('fold change (%s)'%(g1), fontsize=12)
    plt.ylabel('fold change (%s)'%(g2), fontsize=12)

def getAvePoint(ge, gi):
    giAve=[]
    geAve=[]
    geDisc = np.arange(int(min(ge))-1, int(max(ge))+1,0.2)
    for ii in geDisc:
        index_list = [ge.index(val) for val in ge if (val>=ii and val<(ii+0.2))]
        if len(index_list) == 0:
            continue
        else:
            geAve.append(ii)
            giAve.append(np.mean([gi[ix] for ix in index_list]))
            
    return geAve, giAve

### Main Function ###
if __name__ == '__main__':
    spearman_files = glob.glob('result/spearman/all_*.txt')
    
    minVal_dic = {}
    
    for sp_file in spearman_files:
        sp_df = pd.read_table(sp_file, sep='\t', index_col = 0)
        cell = sp_file[len('result/spearman/all_'):-len('_24H.txt')]
        if cell == 'HT29': continue
        if cell == 'PC3': continue
        
        gene = sp_df['ES'].idxmin()
        val = sp_df['ES'].min()
        minVal_dic[cell] = [gene, val]
    
    
    for cell, val in sorted(minVal_dic.items(), key=lambda i: i[1][1], reverse=False):
        print cell, val
    
    plotCnt = 0
    for cell, (gene, val) in sorted(minVal_dic.items(), key=lambda i: i[1][1], reverse=False):
        
        fc_df = pd.read_table('result/foldChange/all_%s_24H.txt'%(cell), sep='\t', index_col = 0)
        cgi_df = pd.read_table('result/cellGrowthInh/%s_24H.txt'%(cell), sep='\t', index_col = 0)
        concat_df = pd.concat([fc_df, cgi_df])
        concat_df = concat_df.ix[:,concat_df.ix[gene]<0.0]
        ge_list = list(concat_df.ix[gene])
        gi_list = list(concat_df.ix['CGI'])
        
        plotCnt += 1
        plt.subplot(2,2,plotCnt)
        plt.plot(ge_list, gi_list, 'o', mfc='none', markersize=8, color = 'dimgray')
        
        if plotCnt == 1:
            # BRD-K49865102|0.04
            plt.plot([-2.11879508195765],[85.8035707425], 'o', mfc='cadetblue', markersize=8, color = 'dimgray')
            
        if plotCnt == 2:
            # BRD-K12184916|0.37
            plt.plot([-2.310290527],[97.8930694479], 'o', mfc='cadetblue', markersize=8, color = 'dimgray')
            
        plt.title('%s\n%s_%s (%s)'%('pValMean', cell, gene, round(val,3)))
        
        plt.subplot(2,2,plotCnt + 2)
        geAve_list, giAve_list = getAvePoint(ge_list, gi_list)
        plt.plot(geAve_list, giAve_list, '-o', color = 'black')
        
        if plotCnt == 2:
            break
    
    plt.show()    
    


        