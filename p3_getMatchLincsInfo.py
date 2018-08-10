'''
Created on Feb 24, 2017

@author: jmjung
'''
import pandas as pd
from commonFunc import getCols

def getLnsFiltered(cellHr_dic):
    for (cell, hr), exp_list in cellHr_dic.items():
        print cell, hr
        col_dList = getCols('../rawData/Broad_LINCS_Level3_INF_mlr12k_n113012x22268_2015-12-31.gct',3,['id','pr_gene_id','pr_gene_symbol'] + exp_list)
    
        lnsFil_file = open('preProcessResult/LINCS/%s_%s_LINCS.txt'%(cell, hr),'w+')
        for sList in col_dList:
            lnsFil_file.write('\t'.join(sList)+'\n')
        
        lnsFil_file.close()

def getMatchExpID():
    ## get matched cell, hr, comp info in LINCS
    
    matCellHrComp_list = []
    matCompCell_df = pd.read_table('preProcessResult/matCompCell_CGPLINCS.txt', sep='\t',index_col = 0)
    for ii in range(len(list(matCompCell_df.index))):
        (cell, hr) = matCompCell_df.iloc[ii,1].split('|')
        compID_list = matCompCell_df.iloc[ii,5].split('|')
        
        for compID in compID_list:
            matCellHrComp_list.append((cell, hr, compID))
        matCellHrComp_list.append((cell, hr, 'DMSO'))
   
    ## get LINCS experiment id for matched cell, hr, comp
    lns_file = open('../rawData/Broad_LINCS_Level3_INF_mlr12k_n113012x22268_2015-12-31.gct')
        
    lineCnt = 0
    data_dList = []
    for line in lns_file.xreadlines():
        lineCnt += 1
        if lineCnt <= 2: continue
        if lineCnt == 16: break
            
        line_list = line.replace('\n','').replace('\r','').split('\t')
        data_dList.append(line_list)
        
    lns_file.close()
       
    
    lns_df = pd.DataFrame(data_dList[1:],columns=data_dList[0]).set_index('id')
    
    cellHr_dic = {}
    for col in list(lns_df):
        cell = lns_df.ix['cl_center_specific_id', col]
        hr = lns_df.ix['sm_time', col]
        comp = lns_df.ix['sm_lincs_id', col]
        
        if comp.upper() == 'DMSO':
            comp = 'DMSO'
           
        if (cell, hr, comp) in matCellHrComp_list:
            if not cellHr_dic.has_key((cell, hr)):
                cellHr_dic[(cell,hr)] = []
            cellHr_dic[(cell,hr)].append(col)
            
    return cellHr_dic

if __name__ == '__main__':
    cellHr_dic = getMatchExpID()
    print 'LINCS file filtering'
    getLnsFiltered(cellHr_dic)
    