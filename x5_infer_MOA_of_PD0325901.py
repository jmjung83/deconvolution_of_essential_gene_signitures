

import pandas as pd
import scipy.stats as sci
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import glob
from decimal import Decimal
import networkx as nx

def makePPIandTransfacFile():
    out_file = open('../rawData/PPIandTRANSFAC.txt','w+')
    print 'get ppi network'
    sym2entz_file = open("../rawData/sym2entz.txt")
     
    entz2sym = {}
    header=True
    for line in sym2entz_file.xreadlines():
        if header: header=False; continue
         
        sym, entz = line.replace('\n','').replace('\r','').split('\t')
         
        if entz != '':
            entz2sym[entz] = sym
     
    sym2entz_file.close()
     
    ####
    read_file = open("../rawData/BIOGRID-ORGANISM-Homo_sapiens-3.4.148.mitab.txt")
     
    for line in read_file.xreadlines():
        if line.startswith('#'): continue
        line_list = line.replace('\n','').split('\t')
 
        if line_list[11] not in ['psi-mi:"MI:0915"(physical association)','psi-mi:"MI:0407"(direct interaction)']: continue
         
        g1 = line_list[0].split('locuslink:')[1]
        g2 = line_list[1].split('locuslink:')[1]
         
        if g1 not in entz2sym.keys() or g2 not in entz2sym.keys():
            continue
         
        out_file.write(entz2sym[g1]+'\t'+entz2sym[g2]+'\tPP\n')
        out_file.write(entz2sym[g2]+'\t'+entz2sym[g1]+'\tPP\n')
     
    read_file.close()    
    print 'get transfac network'
    
    
    read_file = open("../rawData/TRANSFAC.txt")
    
    for line in read_file.xreadlines():
        line_list = line.replace('\n','').replace('\r','').split('\t')
        if len(line_list) != 2:
            continue
        
        [tf, tgs] = line_list
        
        for tg in tgs.split('|'):
            out_file.write(tf+'\t'+tg+'\tPG\n')
        
    read_file.close() 
    out_file.close()


def getNetInfo():
    net = nx.DiGraph()
    in_file = open('../rawData/PPIandTRANSFAC.txt')
    
    for line in in_file.xreadlines():
        [stt, end, info] = line.replace('\n','').split('\t')
        net.add_edge(stt,end,weight=info)
        
    in_file.close()
    return net


if __name__ == '__main__':
    #makePPIandTransfacFile()
    net = getNetInfo()
    stt_list = [sym for sym in net.nodes() if sym.startswith('MAP2K')]
    for stt in stt_list:
        if not nx.has_path(net, stt, 'EIF4E'): continue
        
        shPaths = list(nx.all_shortest_paths(net, stt, 'EIF4E'))
        for shPath in shPaths:
            info_list = []
            for ii in range(len(shPath)-1):
                info_list.append(net[shPath[ii]][shPath[ii+1]]['weight'])
            
            if info_list == ['PP','PP','PG']:
                print '##',shPath 
                    
        
