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
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    seq = ''
    for line in open(sys.argv[1]):
        if line[0] == '>':
            if seq:
                print locus
                print len(seq)
                print seq
                seq = ''
            locus = (line.split(' ')[0])[1:]
        else:
            seq += line.rstrip('*\n')
        #--------------------------------------------
    if seq:
        print locus
        print len(seq)
        print seq
if __name__ == '__main__':
    main()

