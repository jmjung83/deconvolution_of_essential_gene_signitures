import glob
import pandas as pd
import scipy.stats as sci
import statsmodels.sandbox.stats.multicomp
import math
import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import random

PERT_NUM=99999

def get_empPVal(EG_score, perm_values):
    rank = float(sum(x<EG_score for x in perm_values)+1)/float(PERT_NUM+1)
    return rank

def get_EGscore(fc_list, cgi_list):
    cor, p = sci.spearmanr(fc_list, cgi_list)
    EG_score = cor * np.mean(cgi_list)
    
    return EG_score
    
def get_perm_values(negFC, fc_pool, cgi_pool):
    perm_values=[]
    for ii in range(PERT_NUM):
        perm_fc_list  = random.sample(fc_pool, negFC)
        perm_cgi_list = random.sample(cgi_pool,negFC)
        perm_values.append(get_EGscore(perm_fc_list, perm_cgi_list))
    
    return perm_values

cell_num = sys.argv[1]
gr_num   = sys.argv[2]

cell_dic={0:'A375',
          1:'A549',
          2:'BT20',
          3:'HS578T',
          4:'HT29',
          5:'LNCAP',
          6:'MCF7',
          7:'MDAMB231',
          8:'PC3'}

gr_dic={0:'all',
        1:'lmk'}

if __name__ == '__main__':
    cell = cell_dic[int(cell_num)]
    gr = gr_dic[int(gr_num)]
    perm_dic={}
    
    fc_file='./result/foldChange/{}_{}_24H.txt'.format(gr,cell)
    cgi_file = './result/cellGrowthInh/{}_24H.txt'.format(cell)
    
    print(fc_file, "__start__")
    start_time1 = time.time()

    fc_df  = pd.read_table(fc_file, sep='\t', index_col = 0)
    cgi_df = pd.read_table(cgi_file, sep='\t', index_col = 0)

    fc_pool=list(np.array(fc_df.values).flat)
    fc_pool=[val for val in fc_pool if val<0.0]
    cgi_pool=list(np.array(cgi_df.values).flat)
    
    spearman_dic = {}
    
    concat_df = pd.concat([fc_df, cgi_df])
    for gene in list(fc_df.index):
        fc_list  = list(concat_df.loc[gene ,concat_df.loc[gene]<0.0])
        cgi_list = list(concat_df.loc['CGI',concat_df.loc[gene]<0.0])
        
        negFC = len(fc_list)
        if  negFC < 3:
            continue
        
        EG_score = get_EGscore(fc_list, cgi_list)
        
        if negFC not in perm_dic.keys():
            start_time2 = time.time()
            perm_dic[negFC] = get_perm_values(negFC, fc_pool, cgi_pool)
            print(negFC, time.time()-start_time2, "seconds")
        
        emp_Pval = get_empPVal(EG_score, perm_dic[negFC])
        
        spearman_dic[gene] = [len(list(concat_df)), negFC, EG_score, emp_Pval]

    spearman_df = pd.DataFrame.from_dict(spearman_dic, orient='index')
    spearman_df.columns = ['#_of_exp','#_of_exp(fc<0)','EG_score','empirical_pVal']
    
    pd.DataFrame.to_csv(spearman_df, path_or_buf=fc_file.replace('foldChange','spearman'),sep='\t')
    print(fc_file, "__completed__", time.time() - start_time1, "seconds")
        
        