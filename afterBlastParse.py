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
    print >>sys.stderr, ""
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    count = 0 
    key = 0
    for line in open(sys.argv[1]):
        if line[0] == '>':
            key = 1
            savedline = line
        else:
            if key:
                print savedline
                key = 0
                count += 1
            print line
    
    print >>sys.stderr, count
if __name__ == '__main__':
    main()

