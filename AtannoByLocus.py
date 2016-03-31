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

def readAnno(annofile, head=1):
    aDict = {}
    for line in open(annofile):
        if head:
            head -= 1
            continue
        #--------------------------
        key, value = line.strip().split('\t',1)
        if key not in aDict:
            aDict[key] = value
        else:
            print >>sys.stderr, "duplicated key %s" % key
        #------------------------------
    #---------------End reading-------------
    return aDict
#--------------------------------------------

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s locus anno' % sys.argv[0]
        sys.exit(0)
    #--------------------------------------------
    locusFile = sys.argv[1]
    annofile = sys.argv[2]
    aDict = readAnno(annofile)
    for line in open(locusFile):
        locus = line.strip()
        print "%s\t%s" % (locus, aDict[locus])
#-----------------------------------------------------    
if __name__ == '__main__':
    main()

