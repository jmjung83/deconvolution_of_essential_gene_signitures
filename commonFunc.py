'''
Created on Mar 3, 2017

@author: jmjung
'''
import time
import sys

def measure(init_time):
    after_time = time.time()
    dif_time = after_time-init_time
    hour=int(dif_time/3660)
    mins=int((dif_time-hour*3660)/60)
    sec=dif_time-hour*3660-mins*60
    print 'Processing Time :' + str(hour) + 'hour' + str(mins) + 'min' + str(sec) + 'sec'
    
def getCols(fileName, colNameLine, colName_list):
    in_file = open(fileName)
    
    col_dList = []
    
    lineCnt = 0
    data_flg = False
    for line in in_file.xreadlines():
        lineCnt += 1
        #if lineCnt%10000==0:print lineCnt
        line_list = line.replace('\n','').split('\t')
        
        if lineCnt == colNameLine:
            colIndex_list = getColIndex(line_list, colName_list)
            col_dList.append([line_list[i] for i in colIndex_list])
            data_flg = True
            continue
        
        if data_flg:
            col_dList.append([line_list[i] for i in colIndex_list])
        
    in_file.close()
    
    return col_dList

def getColIndex(line_list, colName_list):
    colIndex_list = []
    
    for colName in colName_list:
        try:
            colIndex_list.append(line_list.index(colName))
        except:
            print colName
            raise ValueError ('No exist column names')
            
    return colIndex_list

def getRows(fileName, rowNameLine, rowName_list):
    in_file = open(fileName)
    row_dList = []
    
    rowNameFlag = False
    for line in in_file.xreadlines():
        line_list = line.replace('\n','').split('\t')
        
        if line_list[rowNameLine-1] in rowName_list:
            row_dList.append(line_list[rowNameLine-1:])
            
        if len(row_dList) == len(rowName_list):
            rowNameFlag = True
            break             
    
    if not rowNameFlag:
        raise ValueError ('No exist row names')
        
    in_file.close()
    
    return row_dList

def transpose(dList):
    return map(list, zip(*dList))

def printDic(dic, printCnt):
    cnt = 0
    for key, val in dic.items():
        cnt += 1
        print key, val
        
        if printCnt == cnt: break

def printDList(dList, printCnt):
    cnt = 0
    for sList in dList:
        cnt += 1
        print sList
        
        if printCnt == cnt: break
