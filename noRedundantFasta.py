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
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------
    repDict = readFasta(sys.argv[1])
    #----------------------------------
    tmpSet = set()
    for key, value in repDict.items():
        if value not in tmpSet:
            tmpSet.add(value)
            print '>%s\n%s' % (key, value)
    #------------------------------------------
if __name__ == '__main__':
    main()

