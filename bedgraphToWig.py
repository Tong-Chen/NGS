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
    if len(sys.argv) != 4:
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s bedgraph[ \
usually get using bigwigtowig operate on bigwig file] \
step span' % sys.argv[0]
        sys.exit(0)
    #-------------------------------
    """
    bdg: zero-based,  half-open.
    wig: 1-based, no strand specific
    bed: zero-based,   half-open.
    """
    step = sys.argv[2]
    span = sys.argv[3]
    #-------------------------------------
    for line in open(sys.argv[1]):
        if line[0] == '#':
            lineL = line.split()
            chr, start_end = lineL[2].split(':')
            start = start_end.split('-')[0]
            print "fixedStep chrom=%s start=%s step=%s span=%s" \
                % (chr, start, step, span)
        #---------------------------
        else:
            lineL = line.split()
            if start == lineL[1]:
                print lineL[-1]
                start = lineL[2]
            else:
                start = lineL[1]
                print "fixedStep chrom=%s start=%s step=%s span=%s" \
                    % (lineL[0], start, step, span)
                print lineL[-1]
                start = lineL[2]
        #----------------------------------------------------------------
    #---------------------------------
#----------------------------------------------------------
if __name__ == '__main__':
    main()

