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
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s bdgfile' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------
    """
    bdg: zero-based,  half-open.
    wig: 1-based
    bed: zero-based,   half-open.
    """
    chr = ''
    for line in open(sys.argv[1]):
        lineL = line.strip().split()
        assert(len(lineL) == 4)
        if chr:
            if chr != lineL[0]:
                chr = lineL[0]
                print 'variableStep chrom=%s' % (chr)
        else:
            chr = lineL[0]
            print 'variableStep chrom=%s' % (chr)
        #----------------------------------------
        start = int(lineL[1]) + 1
        end = int(lineL[2]) + 1
        assert(start < end)
        value = lineL[3]
        for i in range(start, end):
            print "%d\t%s" % (i, value)
if __name__ == '__main__':
    main()

