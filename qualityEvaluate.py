#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
import sys
#sys.path.append("/home/ct/pylib")
'''
This file is used to do quality evaluation of the rearranged result.
It draws the picture of the number of locus about the same length 
repeat and different repeat times.

Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
==============================================================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=============================================================================================

if len(sys.argv) != 3:
    print 'Please use %s filename(rearrange) pathToOutput/' % sys.argv[0]
    sys.exit(1)

DEBUG = 0
if DEBUG:
    sum = 0

try:
    fh = open(sys.argv[1])
except IOError, e:
    print 'can not open file %s' % sys.argv[0]
else:
    aDict = {}
    aArray = []
    for line in fh:
        if line.startswith('>'):
            locus = line[1:-1]
            aDict[locus] = {}
            aArray.append(locus)
        else:
            firstSymbol = line.find(':')
            #secondSymbol = line.find(':')
            if firstSymbol != -1:
                length = firstSymbol
                #patched 20110910, before count = line.count(':'), but
                #after filter, there is no '#' for the last rep. 
                count = line.count(':')
                #to decide whether this length of this protein is saved
                if length in aDict[locus]:
                    #to dicide whether this repeat time of this length sequence saved
                    aDict[locus][length].add(count)
                    #to dicide whether this repeat time of this length sequence saved
                else: #this length of this protein is not saved
                    aDict[locus][length] = set([count])
            else:
                print 'wrong'
                sys.exit(1)
    #end reading
    fh.close()
    #-------------------------------------------prase----------------------------------------
    bDict = {}
    for locus, nestDict in aDict.items():
        for length, countArray in nestDict.items():
            for count in countArray:
                if length in bDict:
                    if count in bDict[length]:
                        bDict[length][count].add(locus)
                    else:
                        bDict[length][count] = set([locus])
                else:
                    bDict[length] = {count:set([locus])}
        #end nested dict search for
    #end original dict search for
    #-------------------------length and repeat-time(cnt) distribute-------------------------
    lenDict = {}
    cntDict = {}
    for length, nestDict in bDict.items():
        for count, locusSet in nestDict.items():
            for locus in locusSet:
                if length in lenDict:
                    lenDict[length].add(locus)
                else:
                    lenDict[length] = set([locus])
                if count in cntDict:
                    cntDict[count].add(locus)
                else:
                    cntDict[count] = set([locus])
    #-------------------------------------------output small files-----------------------------
    path = sys.argv[2]
#    for length, nestDict in bDict.items():
#        filename = path + 'len.'+str(length)
#        fh = open(filename, 'w')
#        countArray = nestDict.keys()
#        countArray.sort()
#        for count in countArray:
#            print >>fh, count, len(nestDict[count])
#        fh.close()
    #-------------------------------------------output distribute files--------------------------
    distriLen = path + sys.argv[1] + '.Distri.Length'
    distriCnt = path + sys.argv[1] + '.Distri.Reptime'
    fh = open(distriLen, 'w')
    lenArray = lenDict.keys()
    lenArray.sort()
    for length in lenArray:
        print >>fh, length, len(lenDict[length])
    fh.close()

    fh = open(distriCnt, 'w')
    cntArray = cntDict.keys()
    cntArray.sort()
    for cnt in cntArray:
        print >>fh, cnt, len(cntDict[cnt])
    fh.close()
    

if __name__ == '__main__':
    pass

