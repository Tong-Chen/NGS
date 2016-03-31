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
    print >>sys.stderr, "Print the result to screen"
    print >>sys.stderr, "This substitue the names with AGI in\
downloaded files called aracyc_pathways.20110406"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s alias pathway' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------------------
    symbolL = []
    adict = {}
    head = 1
    for line in open(sys.argv[1]):
        if head:
            head -= 1
            continue
        agi, symbol, other = line.split('\t')
        symbol = symbol.upper()
        agi = agi.upper()
        if symbol in adict:
            if adict[symbol] != agi and \
                adict[symbol].find(symbol) == -1 and \
                agi.find(symbol) == -1:
                agi += '\n' + adict[symbol]
                symbolL.append(symbol)
                print >>sys.stderr, line
                print >>sys.stderr, adict[symbol]
                #return 1
        adict[symbol] = agi
    #----------------------------------------------------
    head = 1
    tr = 0
    for line in open(sys.argv[2]):
        if head:
            print line,
            head -= 1 
            continue
        lineL = line.rstrip().rsplit('\t', 1)
        symbol = lineL[1].upper()
        if symbol != 'UNKNOWN':
            if symbol in adict:
                if symbol in symbolL:
                    print >>sys.stderr, '**', symbol, adict[symbol]
                print '\t'.join((lineL[0], adict[symbol]))
                tr += 1
            else:
                print line,
    #-----------------------------------
    print >>sys.stderr, tr
if __name__ == '__main__':
    main()

