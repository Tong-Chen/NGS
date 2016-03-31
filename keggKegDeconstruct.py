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
    print >>sys.stderr, "This is used to deconstruct the result keg \
files. It outputs many subfiles named by their related level."
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------------------
    file = sys.argv[1]
    for line in open(file):
        if line.startswith('#ENTRY'):
            ko = line.strip().split()[1]
        elif line.startswith('#DEFINITION'):
            define = '_'.join(line.split()[1:-1])
        

if __name__ == '__main__':
    main()

