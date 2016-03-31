#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Transfer locus between tair symbol species.

Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
#from __future__ import division, with_statement
import sys
import re

def main():
    if len(sys.argv) != 3:
        print 'Using python %s alltfgo tair' % sys.argv[0]
        sys.exit(0)
    tair = []
    for line in open(sys.argv[2]):
        tair.append(line[:-1])
    
    symbol = re.compile('symbol:(.+?) ')
    tairc = re.compile('TAIR:(.+?) ')
    aDict = {}
    for line in open(sys.argv[1]):
        if line.startswith('>') and line.find('TAIR:') != -1:
            symbolL = symbol.search(line).group(1)
            tairL = tairc.search(line).group(1).upper()
            aDict[tairL] = symbolL
    
    #--------------get the tair with repetitions-----
    for locus in tair:
        print locus, aDict[locus]

if __name__ == '__main__':
    main()

