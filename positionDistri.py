#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
import sys
from random import randint
from ctIO import readRep, readSeq

def main():
    print >>sys.stderr, "position distribution, divide protein\
into three equal length segments, N-, mid-, C- terminal. Detect by\
midpoint of repetitions. If located in boundary, random choose"
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s seq.for.c rep' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------------------------------
    repDict = {}
    seqDict = {}
    readSeq(sys.argv[1], seqDict)
    readRep(sys.argv[2], repDict)

    for key, valueL in repDict.items():
        length = len(seqDict[key])
        first = length / 3
        second = first * 2
        #print key, first, second, length 
        for valueD in valueL:
            for keys in valueD.keys():
                #print keys,
                mid = sum(keys) / 2
                if mid < first:
                    print -1
                elif mid == first:
                    print -1 if randint(0,1) else 0
                elif mid < second:
                    print 0
                elif mid == second:
                    print 0 if randint(0,1) else 1
                else:
                    print 1
    #----------------------------------------------------

if __name__ == '__main__':
    main()

