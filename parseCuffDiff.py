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

def main():
    if len(sys.argv) != 2:
        print >>sys.stderr, "Transform cuffdiff output."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #----------------------------------------
    head = 1
    aDict = {}
    geneS = set()
    for line in open(sys.argv[1]):
        if head:
            head -= 1
            continue
        #--------------------------------
        lineL = line.strip().split()
        gene = lineL[0]
        key1 = lineL[4]
        key2 = lineL[5]
        if key1 not in aDict:
            aDict[key1] = {}
        if gene not in aDict[key1]:
            aDict[key1][gene] = lineL[7]
            geneS.add(gene)
        if key2 not in aDict:
            aDict[key2] = {}
        if gene not in aDict[key2]:
            aDict[key2][gene] = lineL[8]
    #------------------------------------------------------
    #--------------------------------------------
    allKeyL = aDict.keys()
    allKeyL.sort()
    geneL = list(geneS)
    geneL.sort()
    print "gene\t%s" % ('\t'.join(allKeyL))
    for gene in geneL:
        tmpL = [gene]
        for sample in allKeyL:
            tmpL.append(aDict[sample][gene])
        print "\t".join(tmpL)
#-------------------------------------------------------
if __name__ == '__main__':
    main()

