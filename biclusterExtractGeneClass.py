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
import os

def readMatrix(file, head=1):
    '''
Input file:
    One head or not.
    Two or three columns seperated. Use <datacol> to indict the columns
    containing the data. It must be the last column. So <datacol> is
    always -1.

    #In the data column, <lenD> is used to set the width of each datum.
Outout result:
    aDict = {dataStr:locus}
    '''
    aDict = {}
    for line in open(file):
        if head:
            head -= 1
            continue
        #----------------
        locus, data = line.strip().rsplit('\t', 1)
        #assert lenD > 0, "Too little <lenD>, no data?"
        #assert isinstance(lenD, int), "Not integer number lenD."
        #dataList = [data[i:i+lenD] for i in range(0, len(data), lenD)]
        aDict[data] = locus
    #----------------------------------
    return aDict
#-------------------------------------------------------------------------

def extractLocus(clufile, mDict):
    fh = open(clufile+'.WithLocus', 'w')
    for line in open(clufile): 
        if len(line) < 4: #1 2\n
            break
        line = line.rstrip()
        line = line.replace('-', '')
        if line in mDict:
            print >>fh, "%s\t%s" % (mDict[line], line)
        else:
            print >>sys.stderr, " Wrong line ", line
            sys.exit(1)
    #----------------------------------------------------
    fh.close()
#-------------------------------------------------------------------------

def main():
    print >>sys.stderr, "This is used to integrate gene locus into \
the clustered data. Make sure the length of each column data in \
cluster is the same as in matrix. Another point needs to pay heed to \
is that make sure there is no coulmn deletion or transposition \
occured."
    print >>sys.stderr, "Print the result to files"
    if len(sys.argv) != 5:
        print >>sys.stderr, 'Using python %s matrixname \
matrixHead(0 for no head info, other \
positive integer represents correlate number of heads) \
clusternameLabel(the exact phrase which can represent and can only \
represent one class cluater files. I often use the succeeding parts \
without preceding numbers) clusterpath/' % (sys.argv[0]) 
#diffLen (clusterLenD-matrixLenD; lenD a number of data bits, 5 for %sd, 4 for %sd, 5-4)' % (sys.argv[0], 's', 's')
        sys.exit(0)
    #--------------------------------
    file = sys.argv[1]
    head = int(sys.argv[2])
    cluster = sys.argv[3]
    path = sys.argv[4]

    mDict = readMatrix(file, head)
    #-------------------------------------------
    for clufile in os.listdir(path):
        if clufile.endswith(cluster):
            extractLocus(path+clufile, mDict)
if __name__ == '__main__':
    main()

