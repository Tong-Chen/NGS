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

'''
Input file
Seq 1   2   3   4   5
ACGT    0   1   0   0
AATT    3   0   0   0
TTAA    3   0   0   6
ACTG    0   0   0   1
CAGT    0   1   0   9


output file
Seq 1   2   3   4   5
ACGT    0   1   0   0
AATT    6   0   0   6
ACTG    0   1   0   10


'''



import sys
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"


def rc(seq):
    dict = {'A':'T', 'G':'C', 'T':'A','C':'G','a':'t','g':'c','t':'a','c':'g'}
    seqL = list(seq)
    seqL.reverse()
    tmpL = []
    return ''.join([dict[i] for i in seqL])
#--------------------------------------------------


def main():
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, "This is used for processing motif \
analysis. For usually positive strand and negative strand has no \
difference. This is designed to combine motif and its complementary \
one together."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename header[default 1]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    if lensysargv > 2:
        header = int(sys.argv[2])
    else:
        header = 1
    #-------------------------------------
    col = 0
    aDict = {}
    for line in open(sys.argv[1]):
        if header:
            print line,
            header -= 1
            continue
        #-----------------------------
        lineL = line.split()
        aDict[lineL[col]] = [i for i in lineL[1:]]
    #----------------------------------------------------
    
    #---------------------------------------------------
    keyL = aDict.keys()
    for seq in keyL:
        if seq not in aDict:
            continue
        seq_rc = rc(seq)
        if seq_rc != seq and seq_rc in aDict:
            lenL = len(aDict[seq])
            for i in range(lenL):
                aDict[seq][i] = str(float(aDict[seq][i]) + \
                    float(aDict[seq_rc][i]))
            #-----------------------------------
            aDict.pop(seq_rc)
    #-----------------o--------------------------------------
    #-------new dict and new keys-----
    keyL = aDict.keys()
    keyL.sort()
    for key in keyL:
        print "%s\t%s" % (key, '\t'.join(aDict[key]))
#------------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


