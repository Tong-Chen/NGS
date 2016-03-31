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
from ctIO import readseq, readAnnoNew, readRep
def main():
    print >>sys.stderr, "Paste the mother sequene and related\
 repetition together"
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) < 3:
        print >>sys.stderr, 'Using python %s forc repResult \
[anno] [locus] ' % sys.argv[0]
        sys.exit(0)
    tair = \
        "http://www.arabidopsis.org/servlets/TairObject?name=&type=locus"
    seqDict = {}
    readseq(sys.argv[1], seqDict)
    #--------------------------------------------
    isAnno = 0
    if len(sys.argv) > 3:
        if sys.argv[3]:
            isAnno = 1
            annoDict = {}
            readAnnoNew(sys.argv[3], annoDict)
    #--------------------------------------------
    isLoc = 0
    if len(sys.argv) > 4:
        if sys.argv[4]:
            isLoc = 1
            locL = [locus.strip() for locus in open(sys.argv[4])]
    #--------------------------------------------
    repDict = {}
    '''
    repDict = 
    {
    AT1G62760: 
        [
            {(25, 31): 'SSLSPSS', (51,57): 'SSLSPSS'},
            {(52, 60): 'SLSPSSPPP',},    
        ] 
    }
    '''
    readRep(sys.argv[2], repDict)
    #-------------------------------------------
    repDictKeyL = repDict.keys()
    repDictKeyL.sort()
    for locus in repDictKeyL:
        if isLoc and (locus not in locL) and \
            (locus[:-2] not in locL):
            continue
        #------------------------------------
        print '>', locus
        locusSub = "name=" + locus[:-2]
        print tair.replace("name=", locusSub)
        if isAnno and locus in annoDict:
            print annoDict[locus].replace('\\\\', '\n')
        #if locus in seqDict:
        #    print seqDict[locus]
        seq = seqDict[locus] #this substitute the last one for we
        #know it has this key, if not, wrong
        locusRepL = repDict[locus]
        repoutput = []
        for repitemDIct in locusRepL:
            grep = ''
            posKey = repitemDIct.keys()
            posKey.sort()
            for pos in posKey:
                rep = repitemDIct[pos]
                grep += ':'.join([rep, str(pos[0]), str(pos[1])])
                reps = '*'+rep+'*'
                seq = seq.replace(rep, reps)
            repoutput.append(grep)
        print seq
        print '\n'.join(repoutput)

if __name__ == '__main__':
    main()

