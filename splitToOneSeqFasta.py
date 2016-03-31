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
from ctIO import readFasta

def main():
    print >>sys.stderr, "Print the result to files"
    print >>sys.stderr, "Split a multiple sequence fasta file to\
multiple files with one sequence each"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #---------------
    seqDict = readFasta(sys.argv[1])
    for key, value in seqDict.items():
        fh = open(key, 'w')
        print >>fh, '>%s\n%s' % (key, value)
        fh.close()
if __name__ == '__main__':
    main()

