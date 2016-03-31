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
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
def main():
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, "Transfer 3 column position position value \
    file into a matrix"
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename \
missing_val[default NA, for r plot]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    aDict = {}
    key1L = set()
    key2L = set()
    mv = 'NA'
    if lensysargv >2:
        mv = sys.argv[2]
    for line in open(sys.argv[1]):
        lineL = line.split()
        if len(lineL) == 2:
            lineL.append(mv)
        #---------------------------------
        key1, key2, value = lineL
        key1L.add(key1)
        key2L.add(key2)
        if value != mv and float(value) < 0:
            value = str((-1) * float(value))
        if key1 not in aDict:
            aDict[key1] = {}
        if key2 not in aDict[key1]:
            aDict[key1][key2] = value
        else:
            print >>sys.stderr, "Duplicate key2 %s" % key2
            sys.exit(1)
    #----------------------------------------------
    key1L = list(key1L)
    key1L.sort(key=lambda x: int(x))
    key2L = list(key2L)
    key2L.sort(key=lambda x: int(x))
    print "Sample\t%s" % '\t'.join(key2L)
    for key1 in key1L:
        tmpL = [key1]
        for key2 in key2L:
            tmpL.append(aDict[key1][key2])
        print "\t".join(tmpL)
#-------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


