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
from blast import readMsaPsiBatch

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 5:
        print >>sys.stderr, 'Using python %s msapath/ hitpath/ \
identity[70] evalue[0.01]' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------
    msapath = sys.argv[1]
    hitpath = sys.argv[2]
    identity = float(sys.argv[3])
    evalue = flaot(sys.argv[4])
    psiD = readMsaPsiBatch(hitpath, identity, evalue)
if __name__ == '__main__':
    main()

