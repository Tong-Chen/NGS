#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''


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

def main():
    if len(sys.argv) < 3:
        print >>sys.stderr, "Find common between multiple files and output to the screen"
        print 'Using python %s filename1 filename2 ...' % sys.argv[0]
        sys.exit(0)
    #------------------------------------
    aSet = set([line.strip() for line in open(sys.argv[1])])
    newSet = aSet
    for file in sys.argv[2:]:
        newSet = newSet.intersection( \
            set([line.strip() for line in open(file)]))
    #---------------------------------
    print "\n".join(newSet)

if __name__ == '__main__':
    main()

