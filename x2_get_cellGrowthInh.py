'''
Created on Mar 6, 2017

@author: jmjung
'''

import glob
import pandas as pd
import numpy as np
import math

def getCgpData():
    lns2cgpCell = {}
    lns2cgpComp = {}
    
    in_file = open('preProcessResult/matCompCell_CGPLINCS.txt')
    
    header = True
    for line in in_file.xreadlines():
        if header: header=False;continue
        
        uniCell, cgpCell, lnsCellHr, numComp, pcID, cgpID, lnsID = line.replace('\n','').replace('\r','').split('\t')
        lns2cgpCell[lnsCellHr.split('|')[0]] = cgpCell
        
        lnsID_list = lnsID.split('|')
        cgpID_list = cgpID.split('|')
        
        for ii in range(len(lnsID_list)):
            lns2cgpComp[lnsID_list[ii]] = cgpID_list[ii]
        
    in_file.close()
    
    return lns2cgpCell, lns2cgpComp

def getCellGrowthInhibition(beta, dose, ic50):
    natE = math.e
    t1 = (ic50-dose)*beta
    CGI = 100.0/float(1+math.pow(natE,t1))
    
    return CGI

if __name__ == '__main__':
    lns2cgpCell, lns2cgpComp = getCgpData()
    
    CGI_df = pd.read_table('../rawData/CGP_sdata1_filtered.txt', sep='\t', index_col = 0)
    
    all_files = glob.glob('result/foldChange/all_*.txt')
    for in_file in all_files:
        cell, hr = in_file[len('result/foldChange/all_'):-4].split('_')
        lns_df = pd.read_table(in_file, sep='\t', index_col = 0)
        
        CGI_dic = {}
        CGI_list = []
        
        for compDose in list(lns_df):
            comp, dose = compDose.split('|')
            dose = math.log(float(dose))
            cgpCell = lns2cgpCell[cell]
            cgpComp = lns2cgpComp[comp]
            
            beta = float(CGI_df.ix[cgpCell, cgpComp+'_BETA'])
            ic50 = float(CGI_df.ix[cgpCell, cgpComp+'_IC_50'])
            
            CGI = getCellGrowthInhibition(beta, dose, ic50)
            
            CGI_list.append(CGI)
            
            CGI_dic[compDose] = {'CGI': CGI}
         
        pd.DataFrame.to_csv(pd.DataFrame(CGI_dic), path_or_buf='result/cellGrowthInh/%s_%s.txt'%(cell, hr),sep='\t')
        
        if hr == '24H':
            print cell, hr, len(list(lns_df)), len(CGI_list)
        
        
        
    
    