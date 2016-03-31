#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
import sys


def sortFunc(x, index):
    xl = x.split('\t')
    sum = 0
    for i in index:
        sum += float(xl[int(i)])
    return sum
#-----------------------------------------------------------



def main():
    if len(sys.argv) < 4:
        print >>sys.stderr, "Parse the output of heatmapS.sh to \
reorganize the data. The first column in the head line of file \
is unneeded. That is to say the number of columns of head line is \
1 less than that of data line."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename(cluster.final)[-] \
colSort[3-2-1 means third column be the first in new file, first be \
the third(1-based, only data colum, exclude the first name column. \
Default data value will be sorted by new first data column.)] \
rowSort[4-2-3-5-1 means forth cluster rows be the first 1,from \
down-top] \
clusterSort[Default each cluster is sorted by values in first \
sample in new file. Accept 1-2-3 means first cluster(here cluster 1) \
sort by first \
sample, second cluster sort by second sample. Or 1-2.3-4 means second \
custer sort by the sum of sample2 and sample3. Here the sort of this \
parameter should be accordant with rowSort. Elements separated by dash \
should be the same number of cluster. Actually they are \
negative correlation,  the program will treat with it.' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------------
    file = sys.argv[1]
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    colSort = [int(i) for i in sys.argv[2].split('-')]
    rowSort = [int(i) for i in sys.argv[3].split('-')]
    if len(sys.argv) > 4:
        clusterInnerSort=[i.split('.') for i in sys.argv[4].split('-')]
    else:
        clusterInnerSort=''
    #---------------------------
    #print >>sys.stderr, clusterInnerSort
    head = 1
    allLinesD = {}
    keySet = set()
    for line in fh:
        if head:
            line = 'Gene'+'\t'+line
            head -= 1
        lineL = line.strip().split('\t')
        newline = [lineL[0]]
        for new in colSort:
            newline.append(lineL[new])
        newline = '\t'.join(newline)
        key = lineL[-1]
        keySet.add(key)
        if key not in allLinesD:
            allLinesD[key] = [newline]
        else:
            allLinesD[key].append(newline)
    #---------------------------------------------
    if file != '-':
        fh.close()
    keySetL = list(keySet)
    keySetL.remove('cluster')

    #newKeySetL = [float(i) for i in keySetL if i!= 'cluster']
    keySetL.sort(key=lambda x:float(x))
    print allLinesD['cluster'][0]
    if clusterInnerSort == '':
        for i in rowSort:
            key = keySetL[i-1]
            tmpList = allLinesD[key]
            tmpList.sort(key=lambda x:float(x.split('\t')[1]))
            print "\n".join(tmpList)
    else:
        assert len(clusterInnerSort) == len(rowSort)
        j = -1
        for i in rowSort:
            key = keySetL[i-1]
            tmpList = allLinesD[key]
            tmpList.sort(key=lambda x: sortFunc(x, clusterInnerSort[j]))
            j = j - 1
            print "\n".join(tmpList)

#----------------------------------------------------------
if __name__ == '__main__':
    main()



