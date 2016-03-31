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

def goeast(file, num, head = 1):
    '''
    barDict = {functionname: num}

    heatmapDict={locus: [functionName]}
    '''
    barDict = {}
    barLocus = set()
    heatmapDict = {}
    funcList = []
    for line in open(file):
        if head:
            head -= 1
            continue
        assert line.find('GO:') != -1
        lineL = line.split('\t')
        function, tmpNum, locus = lineL[2], int(lineL[3]), lineL[7]
        #-----bar output----------------------
        numLocus = locus.count(':')
        assert numLocus == tmpNum
        if function not in barDict:
            barDict[function] = numLocus
        else:
            print >>sys.stderr, "Duplicate functiona %s" % function
            sys.exit(1)
        #-----bar output----------------------
        #-----heatmap output------------------
        #function AT1G12 AT2G**
        #cell wall 1     0
        funcList.append(function)
        for locusItem in locus.split(' // '):
            barLocus.add(locusItem)
            locusItem = locusItem.strip()
            if locusItem not in heatmapDict:
                heatmapDict[locusItem] = []
            heatmapDict[locusItem].append(function)
        #-----heatmap output------------------

    #---------------------------------------------
    barDict["Unannotate"] = num - len(barLocus)
    return barDict, heatmapDict, funcList

def bingo(file, num):
    barDict = {}
    barLocus = set()
    heatmapDict = {}
    funcList = []
    wait = 1
    for line in open(file):
        if wait and line.startswith('GO-ID'):
            wait = 0
            continue
        if not wait:
            lineL = line.rstrip().split('\t')
            function, locus = lineL[-2], lineL[-1]
            #---------barDict--------------------------------
            numLocus = locus.count('|')+1
            if function not in barDict:
                barDict[function] = numLocus
            else:
                print >>sys.stderr, "Duplicate function %s"\
                        % function
                sys.exit(1)
            #---------barDict--------------------------------

            #---------heatmapDict--------------------------------
            funcList.append(function)
            for locusItem in locus.split('|'):
                barLocus.add(locusItem)
                if locusItem not in heatmapDict:
                    heatmapDict[locusItem] = []
                heatmapDict[locusItem].append(function)

            #---------heatmapDict--------------------------------
        #----------------------End reading------------
    barDict["Unannotate"] = num - len(barLocus)
    return barDict, heatmapDict, funcList


def main():
    print >>sys.stderr, "this used to parse the gene ontology data\
get from GOEAST or BINGO, getting the data of function distribution."
    print >>sys.stderr, "Print the result to file"
    if len(sys.argv) != 5:
        print >>sys.stderr, 'Using python %s go type[GOEAST, BINGO] \
num{# of genes used for analysis} outputtype[all, bar, heatmap]' % sys.argv[0]
        sys.exit(0)
    #---------------------------------------------
    if sys.argv[2] == "GOEAST":
        barDict, heatmapDict, funcList = \
            goeast(sys.argv[1], int(sys.argv[3]))
    elif sys.argv[2] == "BINGO":
        barDict, heatmapDict, funcList = \
            bingo(sys.argv[1], int(sys.argv[3]))
    else:
        print >>sys.stderr, "Wrong value, should be <GOEAST> or\
<BINGO>"
        sys.exit(1)
    #-----------------------------------------------
    par5 = sys.argv[4]
    funcList.sort()
    if par5 == 'all' or par5 == 'bar':
        fh = open(sys.argv[1]+'.barData', 'w')
        for key in funcList:
            print >>fh, '%s\t%s' % (key, barDict[key])
        print >>fh, "%s\t%s" % ("Unannotate", barDict["Unannotate"])
        fh.close()
    if par5 == 'all' or par5 == 'heatmap':
        fh = open(sys.argv[1]+'.heatmapData', 'w')
        keyList = heatmapDict.keys()
        keyList.sort()
        print >>fh, "function\t%s" % '\t'.join(keyList)
        for funcItem in funcList:
            numL = []
            for key in keyList:
                if funcItem in heatmapDict[key]:
                    numL.append('1')
                else:
                    numL.append('0')
            #-----------------------------------
            print >>fh, "%s\t%s" % (funcItem, '\t'.join(numL))
        fh.close()
if __name__ == '__main__':
    main()

