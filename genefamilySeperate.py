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
from ctIO import readColumns
def main():
    print >>sys.stderr, "Divide genes into different functional \
annotation clusters and save in correlated files."
    print >>sys.stderr, "Print the result to multiple files"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s filename outputdir/' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------------
    aDict = readColumns(sys.argv[1])
    for key, valueS in aDict.items():
        key = key.strip().replace(' ','_')
        file = sys.argv[2] + key
        fh = open(file, 'w')
        for item in valueS:
            if item == 'NULL':
                continue
            item.replace(';', '\n')
            print >>fh, item
        fh.close()
        
if __name__ == '__main__':
    main()

