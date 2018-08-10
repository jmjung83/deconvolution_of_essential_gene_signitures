'''
Created on Mar 22, 2017

@author: jmjung
'''
import glob
import pandas as pd
import scipy.stats as sci
import statsmodels.sandbox.stats.multicomp
import math
import matplotlib.pyplot as plt
import numpy as np
import sys

def distFun((x1,y1), (x2,y2), (x3,y3)): # x3,y3 is the point
    px = x2-x1
    py = y2-y1

    something = px*px + py*py

    u =  ((x3 - x1) * px + (y3 - y1) * py) / float(something)

    if u > 1:
        u = 1
    elif u < 0:
        u = 0

    x = x1 + u * px
    y = y1 + u * py

    dx = x - x3
    dy = y - y3

    # Note: If the actual distance does not matter,
    # if you only want to compare what this function
    # returns to other results of this function, you
    # can just return the squared distance instead
    # (i.e. remove the sqrt) to gain a little performance

    dist = math.sqrt(dx*dx + dy*dy)

    return dist
        
def getVariousCorr(gene, concat_df):
        
    sel_df = concat_df.ix[:,concat_df.ix[gene] > 0.0]
    if len(list(sel_df)) < 3:
        return ('NA','NA')
    
    fc_list = list(sel_df.ix[gene])
    cgi_list = list(sel_df.ix['CGI'])
    
    cor, p = sci.spearmanr(fc_list, cgi_list)
        
    #log10P = -math.log10(0.0000000001+p)*np.sign(cor)
    #log10Pmean = log10P * (np.mean(cgi_list)*0.01)
    corMean = cor * np.mean(cgi_list)
    
    return_list = [len(list(sel_df)), corMean]
    
    #print return_list
    return (return_list, p)
        

if __name__ == '__main__':
    fc_files = glob.glob('result/foldChange/*.txt')
    
    pValue_list = []
    for fc_file in fc_files:
        print fc_file
        
        cgi_file = 'result/cellGrowthInh/' + fc_file[len('result/foldChange/all_'):]

        fc_df = pd.read_table(fc_file, sep='\t', index_col = 0)
        cgi_df = pd.read_table(cgi_file, sep='\t', index_col = 0)
        concat_df = pd.concat([fc_df, cgi_df])
         
        spearman_dic = {}
        
        for gene in list(fc_df.index):
            (spearmanCorr, pValue) = getVariousCorr(gene, concat_df)
            if spearmanCorr == 'NA':
                continue
            spearman_dic[gene] = [len(list(concat_df))] + spearmanCorr
            pValue_list.append(pValue)
            
        spearman_df = pd.DataFrame.from_dict(spearman_dic, orient='index')
        spearman_df.columns = ['#_of_exp','#_of_exp(fc<0)','ES']
        
        pd.DataFrame.to_csv(spearman_df, path_or_buf=fc_file.replace('foldChange','spearman_STOP'),sep='\t')
            
    print min(pValue_list)
