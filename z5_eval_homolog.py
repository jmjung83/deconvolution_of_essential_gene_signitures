


'''
Created on Mar 22, 2017

@author: jmjung
'''

import pandas as pd
import scipy.stats as sci
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import glob
from decimal import Decimal
from scipy.stats import rankdata

cell_dic = {'A375_SKIN':'A375',
            'A549_LUNG':'A549',
            'BT20_BREAST':'BT20',
            'HS578T_BREAST':'HS578T',
            'HT29_LARGE_INTESTINE':'HT29',
            'LNCAPCLONEFGC_PROSTATE':'LNCAP',
            'MCF7_BREAST':'MCF7',
            'MDAMB231_BREAST':'MDAMB231',
            'PC3_PROSTATE':'PC3'}

lns2cow = {'HT29':'HT29_LARGE_INTESTINE',
           'MCF7':'MCF7_BREAST',
           'BT20':'BT20_BREAST',
           'A549':'A549_LUNG'}

def getSigGene_shRNA_2014Cowley(cell):
    cowCell = lns2cow[cell]
    cow_df = pd.read_table('rawData/Achilles_QC_v2.4.3.rnai.Gs.txt', sep='\t', skiprows=2, index_col = 0)
    cow_df = cow_df[['Description',cowCell]]
    cowAve_df = cow_df.groupby('Description').mean()
    
    cowAve_df.columns = ['score']
    return cowAve_df

def getSigGene_A549_shRNA_2008Luo():
    riger_df = pd.read_table('rawData/RIGER_pgshRNA_A549.txt', sep='\t', index_col = 0)
    riger_df.columns = ['score']
    return riger_df

def getSigGene_A375_CRISPR_2014Shalem():
    GSsigGene_dic = {}
    
    allGene_file = open('rawData/CRISPR_A375_all_genes_2014Shalem.txt')
    for line in allGene_file.xreadlines():
        gene_list = line.replace('\n','').replace('\r','').split('/')[1:]
        for gene in gene_list:
            GSsigGene_dic[gene] = 'x'
        
    allGene_file.close()

    sigGene_list = []
    sigGene_file = open('rawData/CRISPR_A375_depletion_genes_2014Shalem.txt')
    for line in sigGene_file.xreadlines():
        gene, score = line.replace('\n','').replace('\r','').split('\t')
        if GSsigGene_dic.has_key(gene):
            GSsigGene_dic[gene] = 'o'
        
    sigGene_file.close()
    gs_df = pd.DataFrame.from_dict(GSsigGene_dic, orient='index')
    gs_df.columns=['cat']
    
    return gs_df


def getHomologInfo():
    homolog_dic = {}
    in_file = open('rawData/homologene.data.txt')
    
    for line in in_file.xreadlines():
        (num, spec, entz, sym) = line.replace('\n','').replace('\r','').split('\t')[:4]
        if spec == '9606':
            homolog_dic[num] = [sym]

    in_file.close()
    
    ###
    in_file = open('rawData/homologene.data.txt')
    
    for line in in_file.xreadlines():
        (num, spec, entz, sym) = line.replace('\n','').replace('\r','').split('\t')[:4]
        if spec == '9606': continue
        if num in homolog_dic.keys():
            homolog_dic[num].append(spec)
    
    in_file.close()
    
    homolog2_dic = {}
    for key, val in homolog_dic.items():
        homolog2_dic[val[0]] = len(val)

    homolog2_df = pd.DataFrame.from_dict(homolog2_dic, orient='index')
    homolog2_df.columns=['homo_count']    

    print 'homolog2_df',len(homolog2_df.index)
    print 'homolog2_df ave.', homolog2_df.mean()['homo_count']
    print(homolog2_df.head())
    
    return homolog2_df

def get_significance(GS, x, pval):
    if GS<x:
        #if pval<0.000001: return '# #'
        if pval<0.0001: return '****'
        if pval<0.001: return '***'
        if pval<0.01: return '**'
        if pval<0.1: return '*'
    return ''

def draw_figure(fileName, TH):
    plt.figure()
    
    df = pd.read_table(fileName, sep='\t', index_col = 0)
    
    bar_list = []
    x = -2
    cellTick = []
    color_list = []
    for row in df.index:
        x += 2
        bar_list.append([x, df.ix[row,'all.Ave'],''])
        
        if df.ix[row,'sigGS.Ave'] != -1:
            x += 1
            bar_list.append([x, df.ix[row,'sigGS.Ave'],''])
            x += 1
            bar_list.append([x, df.ix[row,'sigCP.Ave'], get_significance(df.ix[row,'sigGS.Ave'],df.ix[row,'sigCP.Ave'],df.ix[row,'GS2CP.Pval'])])
            x += 1
            bar_list.append([x, df.ix[row,'sigBT.Ave'], get_significance(df.ix[row,'sigGS.Ave'],df.ix[row,'sigBT.Ave'],df.ix[row,'GS2BT.Pval'])])
            
            color_list += ['gray','cadetblue','orange','lightcoral']
            cellTick.append(x-1.5)
        else:
            x += 1
            bar_list.append([x, df.ix[row,'sigCP.Ave'],''])
            color_list += ['gray','orange']
            cellTick.append(x-0.5)
            
    print color_list
    x_list = []
    val_list = []
    pval_list = []
    for (x, val, pval) in bar_list:
        x_list.append(x)
        val_list.append(val)
        if x in [8,13]:
            plt.text(x-0.45, val+0.8, pval, fontsize=14)
        else:
            plt.text(x-0.45, val+0.02, pval, fontsize=14)

        
    plt.bar(x_list, val_list, color=color_list)
    plt.xticks(cellTick, list(df.index), fontsize=14)
    plt.yticks(fontsize=14)

    plt.ylabel('Averaged homolog count', fontsize=12)
    plt.title('Averaged homolog count of essential genes (TH=%s)'%TH, fontsize=12)
    
    al_patch = mpatches.Patch(color='gray', label='All')
    gs_patch = mpatches.Patch(color='cadetblue', label='PREV')
    cp_patch = mpatches.Patch(color='orange', label='PROP')
    bt_patch = mpatches.Patch(color='lightcoral', label='BOTH')
    
    plt.legend(handles=[al_patch, gs_patch, cp_patch, bt_patch], fontsize=12, loc='lower left',bbox_to_anchor=(1, 0.5))
    plt.ylim(0,22)
    plt.show()
    #plt.savefig(fileName[:-3]+'png', dpi=600)

def getGSsigGene(cell):
    if cell=='A375':
        return getSigGene_A375_CRISPR_2014Shalem()
    
    if cell in ['MCF7', 'BT20', 'HT29']:
        return getSigGene_shRNA_2014Cowley(cell)
    
    if cell == 'A549':
        return getSigGene_A549_shRNA_2008Luo()

def get_resultFile(resultFileName, TH):

    rs_file = open(resultFileName,'w+')
    rs_file.write('\t'.join(['cell','TH',
                             '#all','all.Ave',
                             '#sigGS','sigGS.Ave',
                             '#sigCP','sigCP.Ave','GS2CP.Pval',
                             '#sigBT','sigBT.Ave','GS2BT.Pval']) + '\n')

    ## homolog df
    homo_df = getHomologInfo()
    ## all shared gene
    all = set(homo_df.index)
    allnum = len(all)
    allave = np.mean(homo_df.ix[all,'homo_count'])

    for cell in ['A375','MCF7','BT20','A549', 'MDAMB231','PC3','LNCAP']:
        print '##########', cell, TH
        ## compound df
        cp_df = pd.read_table('result/spearman/lmk_%s_24H.txt'%(cell), sep='\t', index_col = 0)

        ## compound essential gene
        sigCP = set(cp_df.loc[cp_df['empirical_pVal']<=0.00001].index)
        sigCP = sigCP&all
        sigCPnum = len(sigCP)
        sigCPave = np.mean(homo_df.ix[sigCP,'homo_count'])

        ## previous tech. df
        if cell in ['A375','MCF7','BT20','A549']:
            gs_df = getGSsigGene(cell)
        else:
            gs_df = cp_df
        
        if cell in ['A375','MCF7','BT20','A549']:
            ## shRNA essential gene
            if cell == 'A375':
                sigGS = set(list(gs_df[gs_df['cat'] == 'o'].index))
            else:
                ## previeus tech. essential gene threshold
                rankTH = int(len(gs_df)*TH)
                sigGS = set(list(gs_df.sort_values('score').index)[:rankTH])
                
            sigGS=sigGS&all
            sigGSnum = len(sigGS)
            sigGSave = np.mean(homo_df.ix[sigGS,'homo_count'])
            GS2CPpval = sci.ttest_ind(list(homo_df.ix[sigGS,'homo_count']), list(homo_df.ix[sigCP,'homo_count']), equal_var=False)[1] 
            
            ## both essential gene
            sigBT = set(sigGS)&set(sigCP)
            sigBTnum = len(sigBT)
            sigBTave = np.mean(homo_df.ix[sigBT,'homo_count'])
            GS2BTpval = sci.ttest_ind(list(homo_df.ix[sigGS,'homo_count']), list(homo_df.ix[sigBT,'homo_count']), equal_var=False)[1]
            
        else:
            [sigGSnum,sigGSave,GS2CPpval,sigBTnum,sigBTave,GS2BTpval] = [-1]*6
        
        ## write a file
        rs_file.write('\t'.join([cell] + [repr(val) for val in [TH,
                                                                allnum,allave,
                                                                sigGSnum,sigGSave,
                                                                sigCPnum,sigCPave,GS2CPpval,
                                                                sigBTnum,sigBTave,GS2BTpval]]) + '\n')

    rs_file.close()
if __name__ == '__main__':
    TH = 0.1
    fileName = 'result/evaluation/eval_for_landmarkgene/homolog_%s.txt'%(str(TH).replace('.',''))
    #get_resultFile(fileName, TH)
    draw_figure(fileName, TH)
            
    
    

        










     

    