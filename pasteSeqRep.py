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
    print >>sys.stderr, "Paste the mother sequene and related\
        repetition together"
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 4:
        print >>sys.stderr, 'Using python %s forc repResult anno' % sys.argv[0]
        sys.exit(0)
    
    seqDict = {}
    i = 0
    for line in open(sys.argv[1]):
        if i % 3 == 0:
            locus = line.strip()
            seqDict[locus] = ''
        else:
            seqDict[locus] += line
        i += 1
    #--------------------------------------------
    annoDict = {}
    head = 1
    for line in open(sys.argv[3]):
        if head:
            head -= 1
        else:
            
            lineL = line[:-1].split('\t')
            locus = lineL[0]
            if locus not in annoDict:
                annoDict[locus] = set()
            annoDict[locus].add(lineL[1])
            annoDict[locus].add(tuple(lineL[2:4]))
            annoDict[locus].add(tuple(lineL[4:6]))
            annoDict[locus].add(tuple(lineL[6:8]))
            annoDict[locus].add(tuple(lineL[8:10]))
    #--------------------------------------------
    for line in open(sys.argv[2]):
        if line[0] == '>':
            locus = line[1:].split(' ')[0]
            locus = locus.strip()
            print '>', locus
            if locus in annoDict:
                print annoDict[locus]
            if locus in seqDict:
                print seqDict[locus]
        else:
            print line,

if __name__ == '__main__':
    main()

