'''
Created on Mar 3, 2017

@author: jmjung
'''
import pandas as pd

if __name__ == '__main__':
    
    ## compound per cell in CGP
    print 'processing CGP'
    rs_file = open('preProcessResult/compPerCell_CGP.txt','w+')
    
    cgp_df = pd.read_table('../rawData/CGP_sdata1_filtered.txt', sep='\t',index_col = 0)
    
    for cell in list(cgp_df.index):
        validCol_list = list(cgp_df.ix[cell].dropna().index)
        validComp_list = [val.replace('_BETA','') for val in validCol_list if 'BETA' in val]
        rs_file.write(cell + '\t' + '|'.join(validComp_list) + '\n') 
    
    rs_file.close()

    ## compound per cell in LINCS
    print 'processing LINCS'
    # get header info only
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
    
    # get cell, hour, compound info
    rs_file = open('preProcessResult/compPerCell_LINCS.txt','w+')
    cellComp_dic = {}
    
    lns_df = pd.DataFrame(data_dList[1:],columns=data_dList[0]).set_index('id')
    for col in list(lns_df):
        cell = lns_df.ix['cl_center_specific_id', col]
        hr = lns_df.ix['sm_time', col]
        comp = lns_df.ix['sm_lincs_id', col]
        if (cell == '-666') or (hr == '-666') or (comp in ['-666','DMSO']) :
            continue
        else:
            if not cellComp_dic.has_key((cell,hr)):
                cellComp_dic[(cell,hr)] = set()
            cellComp_dic[(cell, hr)].add(comp)
     
    for (cell, hr), comp_set in cellComp_dic.items():
        rs_file.write(cell+'|'+hr+'\t'+'|'.join(list(comp_set))+'\n')
    
    rs_file.close()
    