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
    if len(sys.argv) < 2:
        print '''File Format:
mm9_chr1_0  0
mm9_chr10_0 0
mm9_chr10_1000  0
mm9_chr10_10000 0
mm9_chr10_100000    0
mm9_chr10_1000000   0
mm9_chr10_10000000  0.0464179
mm9_chr10_100000000 0.208881
mm9_chr10_100001000 0.23209
'''
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename \
add_extra_row_to_separate_chr[Default no adding, any non-false value \
means add 10 rows with \
value 0]' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------------
    add = 0
    value = '0'
    getLen = 1
    if len(sys.argv) > 2:
        add = sys.argv[2]
    #-----------------------------------------------
    aDict = {}
    for line in open(sys.argv[1]):
        lineL = line.split()
        if getLen:
            lenLine = len(lineL)
            value = '\t'.join(list(value * (lenLine-1)))
            getLen = 0
        keyL = lineL[0].split('_')
        i_key = keyL[-1]
        o_key = '_'.join(keyL[1:-1])
        o_key.replace('chr', '', 1)
        i_key = int(i_key)
        if o_key not in aDict:
            aDict[o_key] = {}
        if i_key in aDict[o_key]:
            print >>sys.stderr, "Duplicate %d in $s" % (i_key, o_key)
            sys.exit(1)
        aDict[o_key][i_key] = line 
    #---------------------------------------------
    o_keyL = aDict.keys()
    o_keyL.sort()
    for o_key in o_keyL:
        tmpD = aDict[o_key]
        i_keyL = tmpD.keys()
        i_keyL.sort()
        for i_key in i_keyL:
            print tmpD[i_key],
        #-----------------------------
        if add:
            for i in range(10):
                print "add\t%s" % value
#---------------------------------------------
if __name__ == '__main__':
    main()

