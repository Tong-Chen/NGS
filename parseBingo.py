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

def readBinggoOutput(binggoOutput):
    '''
    '''
    annoDict = {}
    fh = open(binggoOutput)
    line = fh.readline()
    while not line.startswith("GO-ID"):
        line = fh.readline()
        if len(line) == 0:
            print >>sys.stderr, 'No label GO-ID.'
            sys.exit(1)
    #-----------------------------------------
    while len(line):
        line = fh.readline()
        if len(line) == 0:
            break
        lineL = line.strip().rsplit('\t', 2)
        key = lineL[1]
        locusL = lineL[2].split('|')
        if key not in annoDict:
            annoDict[key] = locusL
        else:
            print >>sys.stderr, "Duplicated locus %s" % key
    #-----------------------------------------------------
    fh.close()
    return annoDict
#---------------------------------------------------------- 
def main():
    print >>sys.stderr, "This extracts binggo annotated proteins from\
 bingo output and save results in each file."
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s bingoInput binggoOutput' % sys.argv[0]
        sys.exit(0)
    #---------------------------------------------------
    bingoInput = sys.argv[1]
    binggoOutput = sys.argv[2]

    bingoInputList = [line.strip() for line in open(bingoInput)]
    bingoInputSet = set(bingoInputList)

    annoDict = readBinggoOutput(binggoOutput)
    binggoOutputList = [locus for valueL in annoDict.values() for locus
        in valueL]
    binggoOutputSet = set(binggoOutputList)
    #-------------------------------------------
    assert binggoOutputSet.issubset(bingoInputSet)
    unannoSet = bingoInputSet.symmetric_difference(binggoOutputSet)
    unannoList = list(unannoSet)
    unannoList.sort()
    fh = open('bingoUnannotated.Iso', 'w')
    print >>fh, '\n'.join(unannoList)
    fh.close()
    #-------------------------------------------
    for key, valueL in annoDict.items():
        #if len(valueL) == 0:
        #    print >>sys.stderr, "Wrong valueL %s" % key
        #    sys.exit(1)
        key = key.replace(' ','_')
        key = key.replace('/','_')
        key = key.replace(',','')
        fh = open(key+'.Iso', 'w')
        valueL.sort()
        print >>fh, '\n'.join(valueL)
        fh.close()
    #-------------------------------------------
if __name__ == '__main__':
    main()

