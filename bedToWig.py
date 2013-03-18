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

def output(chr, adict, col):
    '''
    adict={chr:{pos:value}}
    '''
    print 'variableStep chrom=%s' % chr
    posL = adict[chr].keys()
    posL.sort(key=lambda x:int(x))
    if col == -1:
        for i in posL:
            print "%s\t%d" % (i,  adict[chr][i])
    else:
        for i in posL:
            print "%s\t%f" % (i,  adict[chr][i])
    #----------------------------------------
#-----------------------------

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s bedfile[sort by \
chromosome] innerValue[Default FALSE,which means value is \
1. Accept a number to represent the column[1-based] of \
the innerValue.]' % sys.argv[0]
        sys.exit(0)
    #-------------------------------
    """
    bdg: zero-based,  half-open.
    wig: 1-based, no strand specific
    bed: zero-based,   half-open.
    """
    col = -1
    if len(sys.argv) == 3:
        col = int(sys.argv[2]) - 1
    #-------------------------------------
    chr = ''
    adict = {}
    for line in open(sys.argv[1]):
        lineL = line.split()
        if col == -1:
            value = 1
        else:
            value = float(lineL[col])
        #----------------------------------
        if chr and chr != lineL[0]:
            output(chr, adict, col)
            adict = {}
            chr = ''
        if not chr:
            chr = lineL[0]
            assert( chr not in adict)
            adict[chr] = {}
        #----------------------------
        start = int(lineL[1]) + 1
        end = int(lineL[2]) + 1
        for i in range(start,  end):
            if i not in adict[chr]:
                adict[chr][i] = value
            else:
                adict[chr][i] += value
        #---------------------------
    #---------------------------------
    if chr:
        output(chr, adict, col)
if __name__ == '__main__':
    main()

