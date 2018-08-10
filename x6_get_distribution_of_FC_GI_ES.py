'''
Created on Mar 22, 2017

@author: jmjung
'''
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.kde import gaussian_kde
from numpy import linspace

if __name__ == '__main__':
    
    ## fc
    fc_files = glob.glob('result/foldChange/all_*_24H.txt')
      
    totFC = []
    for fc_file in fc_files:
        if ('HT29' in fc_file) or ('HS578T' in fc_file):
            continue 
        fc_df = pd.read_table(fc_file, sep='\t', index_col = 0)
        for col in list(fc_df):
            totFC += list(fc_df[col])
      
    fc_pdf = gaussian_kde(totFC)
    x=linspace(-3.5, 3.5,200)
    plt.plot(x,fc_pdf(x),color = 'crimson')
    plt.hist(totFC,normed=True, bins=200, color='cadetblue')
    plt.xlim(-3.5, 3.5)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.show()
    
    ## cgi
    cgi_files = glob.glob('result/cellGrowthInh/*_24H.txt')
    
    totCGI = []
    for cgi_file in cgi_files:
        if ('HT29' in cgi_file) or ('HS578T' in cgi_file):
            continue 
        cgi_df = pd.read_table(cgi_file, sep='\t', index_col = 0)
        totCGI += list(cgi_df.ix['CGI'])
        
    cgi_pdf = gaussian_kde(totCGI)
    x=linspace(min(totCGI), max(totCGI),200)
    plt.plot(x,cgi_pdf(x),color = 'crimson')
    plt.hist(totCGI,normed=True, bins=100, color='cadetblue')
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.show()

    ## EG score
    eg_files = glob.glob('result/spearman/all_*_24H.txt')
    
    totEG = []
    for eg_file in eg_files:
        if ('HT29' in eg_files) or ('HS578T' in eg_files):
            continue 
        eg_df = pd.read_table(eg_file, sep='\t', index_col = 0)
        totEG += list(eg_df['EG_score'])
        
    eg_pdf = gaussian_kde(totEG)
    #x=linspace(min(totEG), max(totEG),200)
    x=linspace(-40, 40,200)
    plt.plot(x,eg_pdf(x),color = 'crimson')
    plt.hist(totEG,normed=True, bins=150, color='cadetblue')
    plt.xlim(-40, 40)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.show()