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
    print >>sys.stderr, "Transfer homolgo for c file to multiple\
simple fatsa named after the first locus."
    print >>sys.stderr, "Print the result to several files"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #----------------------------------
    fh = ''
    for line in open(sys.argv[1]):
        if line[0] == '=':
            file = line[1:].split()[0]
            if fh:
                fh.close()
            #--------------------------------
            fh = open(file, 'w')
        #----------------------------
        else:
            print >>fh, line,
    #----close for the last file-------
    fh.close() 
if __name__ == '__main__':
    main()

