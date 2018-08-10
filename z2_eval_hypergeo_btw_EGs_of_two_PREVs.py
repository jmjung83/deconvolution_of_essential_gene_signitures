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

def getSigGene_A375_CRISPR_2014Shalem():
    GSsigGene_dic = {}
    
    allGene_file = open('../rawData/CRISPR_A375_all_genes_2014Shalem.txt')
    for line in allGene_file.xreadlines():
        gene_list = line.replace('\n','').replace('\r','').split('/')[1:]
        for gene in gene_list:
            GSsigGene_dic[gene] = 'x'
        
    allGene_file.close()

    sigGene_list = []
    sigGene_file = open('../rawData/CRISPR_A375_depletion_genes_2014Shalem.txt')
    for line in sigGene_file.xreadlines():
        gene, score = line.replace('\n','').replace('\r','').split('\t')
        if GSsigGene_dic.has_key(gene):
            GSsigGene_dic[gene] = 'o'
        
    sigGene_file.close()
    GSsigGene_df = pd.DataFrame.from_dict(GSsigGene_dic, orient='index')
    GSsigGene_df.columns=['cat']
    
    return GSsigGene_df

def getSigGene_shRNA_2010Kim(cell):
    GSsigGene_dic = {}
    
    in_file = open('rawData/%s_shRNA_2010Kim.txt'%(cell))
    
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg = False; continue;
        
        line = line.replace("\n","").replace("\r","").replace('"','')
        line_list = line.split("\t")
        
        gene = line_list[3]
        GSsigGene_dic[gene] = 'x'
        
        EG_info  = line_list[9].split(";")[0]
        if EG_info == "essential gene":
            GSsigGene_dic[gene] = 'o'
            
    in_file.close()
    
    GSsigGene_df = pd.DataFrame.from_dict(GSsigGene_dic, orient='index')
    GSsigGene_df.columns=['cat']
    
    return GSsigGene_df

def isFloat(string):
    try:
        float(string)
    except:
        return False
    else:
        return True

def getSigGene_PC3_shRNA_2012Ros():
    GSsigGene_dic = {}
    
    in_file = open('../rawData/PC3_AshRNA_metabolite_cellMass.txt')
                        
    header = True
    for line in in_file.xreadlines():
        if header: header = False;continue;
             
        line = line.replace('\n','').replace('\r','')
        line_list = line.split(" ")
             
        if not isFloat(line_list[-1]): continue
        gene = line_list[-8].upper()
        
        if (gene == "POOL") or ("#" in gene) or (gene == "RISC-FREE"): continue
        GSsigGene_dic[gene] = 'x'
        
        val_FM = float(line_list[-2])
        val_LM = float(line_list[-1])
        
        if (float(val_FM) <= 0.6) or (float(val_LM) <= 0.6):
            GSsigGene_dic[gene] = 'o'
         
    in_file.close()
    
    GSsigGene_df = pd.DataFrame.from_dict(GSsigGene_dic, orient='index')
    
    GSsigGene_df.columns=['cat']
    return GSsigGene_df

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


def hyperGeoTest(GSsigGene_df, cell, GStype, hyperGeo_dic,hyperGeoPval_dic, TH):
         
    cpOri_df = pd.read_table('result/spearman/all_%s_24H.txt'%(cell), sep='\t', index_col = 0)

    all = list(set(GSsigGene_df.index)&set(cpOri_df.index))
    gs_df = GSsigGene_df.ix[all]
    cp_df = cpOri_df.ix[all]
    
    if list(gs_df) == ['score']:
        rankTH = int(len(list(gs_df.index))*TH*0.01)
        sigGS = list(gs_df.sort_values('score').index)[:rankTH]
    elif list(gs_df) == ['cat']:
        sigGS = list(gs_df[gs_df['cat'] == 'o'].index)    

    

    for measure in ['cor','corMean','log10P','log10Pmean', 'lnP','lnPmean']:
    
        rankTH = int(len(list(cp_df.index))*TH*0.01)
        sigCP = list(cp_df.sort_values(measure).index)[:rankTH]
        hit = set(sigGS)&set(sigCP)
        p_val = sci.hypergeom.sf(len(hit), len(all), len(sigGS), len(sigCP))
        perc_sigCP = len(sigCP)/float(len(all))
        hyperGeo_dic[(cell, GStype, measure, TH)] = [p_val, len(hit), len(all), len(sigGS), len(sigCP), perc_sigCP]
        
        if (cell, GStype, TH) not in hyperGeoPval_dic.keys():
            hyperGeoPval_dic[(cell, GStype, TH)] = []
        hyperGeoPval_dic[(cell, GStype, TH)].append(p_val)


lns2cow = {'HT29':'HT29_LARGE_INTESTINE',
           'MCF7':'MCF7_BREAST',
           'BT20':'BT20_BREAST',
           'A549':'A549_LUNG'}


def getSigGene_comp(cell):
    cp_df = pd.read_table('result/spearman/all_%s_24H.txt'%(cell), sep='\t', index_col = 0)
    return cp_df[['ES']]

if __name__ == '__main__':

    TH=0.1
    
    for cell in ['A549','MCF7']:
        print '########## %s ############'%cell
        if cell == 'A549':
            ref1_df = getSigGene_shRNA_2014Cowley(cell);   ref1 = 'shRNA_genome_Cowley'
            ref2_df = getSigGene_A549_shRNA_2008Luo();       ref2 = 'shRNA_genome_Luo'
            #comp_df = getSigGene_comp(cell);               comp = 'compound_genome_jung'
            
        if cell == 'MCF7':
            ref1_df = getSigGene_shRNA_2014Cowley(cell);   ref1 = 'shRNA_genome_Cowley'
            ref2_df = getSigGene_shRNA_2010Kim(cell);       ref2 = 'shRNA_kinase_Kim'
            print len(set(ref2_df[ref2_df['cat'] == 'o'].index))
            #comp_df = getSigGene_comp(cell);               comp = 'compound_genome_jung'
        
        print ref1, len(ref1_df.index)
        print ref2, len(ref2_df.index)
        
        all = set(ref1_df.index)&set(ref2_df.index)#&set(comp_df.index)
        rankTH = int(len(all)*TH)
        
        ref1_df = ref1_df.ix[all]
        sigRef1 = set(ref1_df.sort_values('score').index[:rankTH])
        
        ref2_df = ref2_df.ix[all]
        if list(ref2_df) == ['score']:
            sigRef2 = set(ref2_df.sort_values('score').index[:rankTH])
        elif list(ref2_df) == ['cat']:
            sigRef2 = set(ref2_df[ref2_df['cat'] == 'o'].index)  

        #comp_df = comp_df.ix[all]
        #sigComp = set(comp_df.sort_values('ES').index[:rankTH])
        
        print 'all', len(all)
        print 'sigRef1', len(sigRef1)
        print sigRef1
        print 'sigRef2', len(sigRef2)
        print sigRef2
        #print 'sigComp', len(sigComp)
        #print sigComp
        
        hit = sigRef1&sigRef2
        print 'sigRef1 (%s) & sigRef2 (%s)'%(ref1, ref2), len(hit)
        print sci.hypergeom.sf(len(hit), len(all), len(sigRef1), len(sigRef2))
        
#         hit = sigRef1&sigComp
#         print 'sigRef1 (%s) & sigComp (%s)'%(ref1, comp), len(hit)
#         print sci.hypergeom.sf(len(hit), len(all), len(sigRef1), len(sigComp))
#          
#         hit = sigRef2&sigComp
#         print 'sigRef2 (%s) & sigComp (%s)'%(ref2, comp), len(hit)
#         print sci.hypergeom.sf(len(hit), len(all), len(sigRef2), len(sigComp))
#         
#         hit = sigRef1&sigRef2&sigComp
#         print 'sigRef1 (%s) & sigRef2 (%s) & sigComp (%s)'%(ref1, ref2, comp), len(hit)
        



        
