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
from ctIO import readseq, readAnno
def main():
    print >>sys.stderr, "Paste the mother sequene and related\
 repetition together"
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) < 3:
        print >>sys.stderr, 'Using python %s forc repResult \
[anno] [locus] ' % sys.argv[0]
        sys.exit(0)
    
    seqDict = {}
    readseq(sys.argv[1], seqDict)
    #--------------------------------------------
    isAnno = 0
    if sys.argv[3]:
        isAnno = 1
        annoDict = {}
        readAnno(sys.argv[3], annoDict)
    #--------------------------------------------
    isLoc = 0
    if sys.argv[4]:
        isLoc = 1
        locL = [locus.strip() for locus in open(sys.argv[4])]
    #--------------------------------------------
    for line in open(sys.argv[2]):
        if line[0] == '>':
            output = 1
            locus = line[1:].strip()
            if isLoc and \
                (locus not in locL) and \
                (locus[:-2] not in locL):
                output = 0
                continue
            #------------------------------------
            print '>', locus
            if isAnno and locus in annoDict:
                print annoDict[locus]
            if locus in seqDict:
                print seqDict[locus]
        else:
            if output:
                print line,

if __name__ == '__main__':
    main()

