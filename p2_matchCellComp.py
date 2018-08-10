'''
Created on Mar 3, 2017

@author: jmjung
'''
#from prevalentFunc import *
import pandas as pd
from commonFunc import printDic

def getCellCompDic(name):
    cellComp_dic = {}
    in_file = open('preProcessResult/compPerCell_%s.txt'%(name))
     
    for line in in_file.xreadlines():
        cell, comp = line.replace('\n','').replace('\r','').split('\t')
        cellComp_dic[cell] = comp.split('|')
         
    in_file.close()
    return cellComp_dic

def getCgpPCDic():
    cgp2pc = {}
    pc2cgp = {}
    
    cgp_df = pd.read_table('../rawData/CGP_sdata2.txt', sep='\t', index_col = 0)
    for cgpID in list(cgp_df.index):
        pc = str(cgp_df.ix[cgpID, 'Pubchem ID'].replace('"','').strip())
        if pc =="NO MATCH": continue
        cgp2pc[str(cgpID)] = pc
        pc2cgp[pc] = str(cgpID) 
        
    return cgp2pc, pc2cgp

def getLnsPCDic():
    # lsm to pc
    in_file = open('../rawData/LSM2PubchemCID.txt')
    
    lsm2pc = {}
    pc2lsm = {}
    header = True
    for line in in_file:
        if header: header=False; continue
        
        if len(line.replace('\n','').split('\t')) != 6: continue
        
        (lsm, pc, smiles, inchi, inchikey, mass) = line.replace('\n','').split('\t')
        if (pc==''): continue
        
        lsm = lsm.strip(); pc = pc.strip()
        lsm2pc[lsm] = pc
        pc2lsm[pc]=lsm
    
    in_file.close()
    
    # brd to pc
    in_file = open('../rawData/BRD2LSM.txt')
    
    brd2pc = {}
    lsm2brd = {}
    header = True
    for line in in_file:
        if header: header=False; continue
        (lsm, brd) = line.replace('\n','').split('\t')
        lsm = lsm.strip()
        brd = brd.strip()
        if not brd.startswith('BRD'): continue
        
        lsm2brd[lsm] = brd
        if lsm2pc.has_key(lsm):
            brd2pc[brd] = lsm2pc[lsm]
    
    in_file.close()
    
    pc2brd = {}
    for pc, lsm in pc2lsm.items():
        if lsm2brd.has_key(lsm):
            pc2brd[pc] = lsm2brd[lsm]
            
    return brd2pc, pc2brd

def unified(cell):
    cell = cell.lower().replace('-','').replace('.','').strip()
    if cell =='lncapclonefgc': cell = 'lncap'
    return cell

def getUnifiedCell(dic, src):
    uni2cell = {}
    if src == 'CGP':
        for key in dic.keys():
            uniKey = unified(key)
            uni2cell[uniKey] = key
    if src == 'LINCS':
        for key in dic.keys():
            uniKey = unified(key.split('|')[0].replace('.101','').replace('.311',''))
            
            if not uni2cell.has_key(uniKey):
                uni2cell[uniKey] = []
            uni2cell[uniKey].append(key)
    
    return uni2cell

if __name__ == '__main__':
    # get pubchem id of CGP compound
    cgpCellComp_dic = getCellCompDic('CGP')
    cgp2pc, pc2cgp = getCgpPCDic()
    
    for cell, comp_list in cgpCellComp_dic.items():
        newComp_list = []
        for comp in comp_list:
            if cgp2pc.has_key(comp):
                newComp_list.append(cgp2pc[comp])
        cgpCellComp_dic[cell] = newComp_list
    
    printDic(cgpCellComp_dic,5)
    
    # get pubchem id of LINCS compound
    lnsCellComp_dic = getCellCompDic('LINCS')
    brd2pc, pc2brd  = getLnsPCDic()
    
    for cellHr, comp_list in lnsCellComp_dic.items():
        newComp_list = []
        for comp in comp_list:
            if brd2pc.has_key(comp):
                newComp_list.append(brd2pc[comp])
        lnsCellComp_dic[cellHr] = newComp_list
    
    printDic(lnsCellComp_dic,5)
    
    # get matched compound per cell
    uni2cgp = getUnifiedCell(cgpCellComp_dic, 'CGP')
    uni2lns = getUnifiedCell(lnsCellComp_dic, 'LINCS')
    
    rs_file = open('preProcessResult/matCompCell_CGPLINCS.txt','w+')
    rs_file.write('\t'.join(['unified_cell', 'CGP_cell','LINCS_cell', '#_of_compound', 'PC_id', 'CGP_id','LINCS_id'])+'\n')
    
    for uniCell in set(uni2cgp.keys()) & set(uni2lns.keys()):
        cgpCell = uni2cgp[uniCell]
        cgpComp_list = cgpCellComp_dic[cgpCell]
        for lnsCell in uni2lns[uniCell]:
            lnsComp_list = lnsCellComp_dic[lnsCell]
            matchComp = list(set(cgpComp_list)&set(lnsComp_list))
            if len(matchComp):
                rs_file.write('\t'.join([uniCell,
                                         cgpCell,
                                         lnsCell,
                                         repr(len(matchComp)),
                                         '|'.join(matchComp),
                                         '|'.join([pc2cgp[pc] for pc in matchComp]),
                                         '|'.join([pc2brd[pc] for pc in matchComp])]) +'\n')
    rs_file.close()