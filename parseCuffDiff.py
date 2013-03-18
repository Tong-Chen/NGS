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
This transfers Cuffdiff output to normal table format to show the gene
expression from each sample.

gene    Brain    Hand    Foot
XLOC_000001     0       0       3.10751
XLOC_000002     57.6877 68.9117 55.3189
XLOC_000003     16.53   15.9654 12.2609
XLOC_000004     0       0.027853        0
XLOC_000005     5.74841 4.48476 6.52879
XLOC_000006     0       0.250831        0
XLOC_000007     0.0109654       0.00795572      2.57567
XLOC_000008     2499.89 2604.53 7591.29
XLOC_000009     50.5431 35.8262 69.7913
'''



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

