'''
Created on Mar 22, 2017

@author: jmjung
'''

import pandas as pd
import scipy.stats as sci
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import glob
import math
import matplotlib.patches as mpatches

def getSigGene_A375_CRISPR_2014Shalem():
    GSsigGene_dic = {}
    
    allGene_file = open('./rawData/CRISPR_A375_all_genes_2014Shalem.txt')
    for line in allGene_file.xreadlines():
        gene_list = line.replace('\n','').replace('\r','').split('/')[1:]
        for gene in gene_list:
            GSsigGene_dic[gene] = 'x'
        
    allGene_file.close()

    sigGene_list = []
    sigGene_file = open('./rawData/CRISPR_A375_depletion_genes_2014Shalem.txt')
    for line in sigGene_file.xreadlines():
        gene, score = line.replace('\n','').replace('\r','').split('\t')
        if GSsigGene_dic.has_key(gene):
            GSsigGene_dic[gene] = 'o'
        
    sigGene_file.close()
    gs_df = pd.DataFrame.from_dict(GSsigGene_dic, orient='index')
    gs_df.columns=['cat']
    
    return gs_df

def getSigGene_shRNA_2014Cowley(cell):
    cowCell = lns2cow[cell]
    cow_df = pd.read_table('./rawData/Achilles_QC_v2.4.3.rnai.Gs.txt', sep='\t', skiprows=2, index_col = 0)
    cow_df = cow_df[['Description',cowCell]]
    cowAve_df = cow_df.groupby('Description').mean()
    
    cowAve_df.columns = ['score']
    return cowAve_df

def getSigGene_A549_shRNA_2008Luo():
    riger_df = pd.read_table('./rawData/RIGER_pgshRNA_A549.txt', sep='\t', index_col = 0)
    riger_df.columns = ['score']
    return riger_df


def hyperGeoTest(gs_df, cell, GStype, hyperGeo_dic, TH):
    cp_df = pd.read_table('result/spearman/lmk_%s_24H.txt'%(cell), sep='\t', index_col = 0)
    sigCP = set(cp_df.loc[cp_df['empirical_pVal']<=0.00001].index)
    
    if list(gs_df) == ['score']:
        rankTH = int(float(len(gs_df))*TH)
        sigGS = set(gs_df.sort_values('score').index[:rankTH])
    elif list(gs_df) == ['cat']:
        sigGS = set(gs_df[gs_df['cat'] == 'o'].index)  

    all = set(gs_df.index)&set(cp_df.index)
    sigCP = all&sigCP
    sigGS = all&sigGS
    hit = sigGS&sigCP

    p_val = sci.hypergeom.sf(len(hit), len(all), len(sigGS), len(sigCP))
    hyperGeo_dic[(cell, GStype, TH)] = [float('%.3g'%p_val), len(hit), len(all), len(sigGS), len(sigCP)]
    
lns2cow = {'HT29':'HT29_LARGE_INTESTINE',
           'MCF7':'MCF7_BREAST',
           'BT20':'BT20_BREAST',
           'A549':'A549_LUNG'}

def getGS_sigGene(cell):
    if cell=='A375':
        return getSigGene_A375_CRISPR_2014Shalem(), 'pgCRISPR_14Shalem'
    
    if cell in ['MCF7', 'BT20']:
        return getSigGene_shRNA_2014Cowley(cell), 'pgshRNA_14Cowley'
    
    if cell == 'A549':
        return getSigGene_A549_shRNA_2008Luo(), 'pgshRNA_08Luo'

    return gs_df, gsType

def draw_figure(dic):
    # hyperGeo_dic[(cell, GStype, TH)] = [p_val, len(hit), len(all), len(sigGS), len(sigCP)]
    x_list = []
    y_list = []
    
    x_pos=[1,2,3,4]
    hit_list=[]
    pval_list=[]
    for key, val in sorted(hyperGeo_dic.items(), key=lambda i: i[1][0], reverse=False):
        print(key, val)
        x_list.append(key[0])
        pval = -math.log10(val[0])
        if pval >= 6: pval = 6
        y_list.append(pval)
        hit_list.append(val[1])
        pval_list.append(val[0])
         
    plt.figure()
    plt.bar(x_pos, y_list, color=['white'], edgecolor=['black'])
    plt.xticks(x_pos, x_list, fontsize=14)
    plt.yticks([0,1,2,3,4,5,6,6.2], ['0','1','2','3','4','6','> 6'], fontsize=12)
    plt.ylabel('-log(p-value)', fontsize=14)
    
    for x,y, pval, hit in zip(x_pos, y_list, pval_list,hit_list):
        if x in [3,4]:
            plt.text(x-0.3,y+0.6,'P:' + str(pval),fontsize=14)
            plt.text(x-0.3,y+0.2,'S:' + str(hit),fontsize=14)
        else:
            plt.text(x-0.3,y-0.4,'P:' + str(pval),fontsize=14)
            plt.text(x-0.3,y-0.8,'S:' + str(hit),fontsize=14)
        
    plt.show()

if __name__ == '__main__':
    
    rs_file = open('result/evaluation/eval_for_landmarkgene/hypergeo_btw_PREV_and_PROP.txt','w+')
    rs_file.write('\t'.join(['cell','GStype','TH','p_val', 'hit', 'all', 'sigGS', 'sigCP']) + '\n')

    hyperGeo_dic = {}
    
    ## rank based
    TH=0.1
    for cell in ['A375','MCF7','BT20','A549']:
        gs_df, gsType = getGS_sigGene(cell)
        hyperGeoTest(gs_df, cell, gsType, hyperGeo_dic,TH)

    for key, val in sorted(hyperGeo_dic.items(), key=lambda i: i[1][0], reverse=False):
        rs_file.write('\t'.join([str(ii) for ii in list(key)+val]) + '\n')

    rs_file.close()
    print(hyperGeo_dic)
    draw_figure(hyperGeo_dic)
    
        
