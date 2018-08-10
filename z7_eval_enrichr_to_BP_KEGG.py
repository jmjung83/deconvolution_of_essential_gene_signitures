'''
Created on 2018. 7. 25.

@author: jmjung
'''
import gseapy

gSet_dic={'KEGG2016':'KEGG_2016',
          'BP2018'  :'GO_Biological_Process_2018',
          'MF2018'  :'GO_Molecular_Function_2018',
          'WP2016'  :'WikiPathways_2016',
          'RE2016'  :'Reactome_2016'}

sig_gene_list = ['POP1','FTSJ2','SYNCRIP','E2F1','SFPQ','BRCA1','DTYMK','DTL','RFC2','DSN1','UNG','DSCC1','UMPS','RFC3','RFC4','GEMIN2','SNRPD1','NCBP1','DONSON','RFWD3','UCHL5','PAICS','TBCCD1','UBE2N','FEN1','TIMELESS','MIS18A','NXT1','SEPHS1','ERCC6L','CHEK2','CHAF1B','CHAF1A','EXO1','PTGES3','POLR3K','EXOSC2','EXOSC9','EIF4E','RAD51','TMEM48','FANCA','MLF1IP','MYBL2','FANCC','SETD8','EEF1E1','EED','PSMC3IP','RBL1','TAF1A','MYO19','RNASEH2A','POLE3','GINS2','TRA2B','ANAPC10','SKP2','TRIP13','MAK16','HSPA14','LSM6','TMEM109','CTPS','CD3EAP','HEATR3','LMNB2','PKMYT1','BLM','RPP30','CSE1L','POLR2D','RQCD1','COX4NB','PPAT','HNRNPAB','ZWILCH','C17ORF42','MCM6','CDC25A','TIPIN','TUBG1','GINS3','ABCE1','C1ORF135','PRIM1','SNRPA1','TYMS','SRSF2','MCM5','MCM4','SRSF3','MCM3','CLN6','MCM10','DCLRE1B','MCM2','NIP7','POLA2','CDC6','SNRNP40','ASF1B','PUS7L']
print(len(sig_gene_list))
gseapy.enrichr(gene_list=sig_gene_list, gene_sets='KEGG_2016', outdir='enrichR/KEGG2016')
gseapy.enrichr(gene_list=sig_gene_list, gene_sets='GO_Biological_Process_2018', outdir='enrichR/BP2018')
