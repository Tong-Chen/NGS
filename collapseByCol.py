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
    print >>sys.stderr, "This collapses each column by the specified \
column[default the first column, 1-based count]. Duplicated items are \
merged together."
    print >>sys.stderr, "Print the result to screen"
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, 'Using python %s filename column[1-based, \
default 1] header[number of header lines, default 0]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    col = int(sys.argv[2])-1 if lensysargv == 3 else 0
    header = int(sys.argv[3]) if lensysargv == 4 else 0
    aDict = {}
    for line in open(sys.argv[1]):
        if header:
            print line,
            header -= 1
            continue
        #--------context--------------------------
        lineL = line.strip('\n').split('\t')
        key = lineL[col]
        lenLineL = len(lineL)
        if key not in aDict:
            aDict[key] = [set([i]) for i in lineL]
        else:
            for i in range(lenLineL):
                aDict[key][i].add(lineL[i])
        #--------------------------------------------
    #-----------------------------------------------
    #-------output------------------------------------
    keyL = aDict.keys()
    keyL.sort()
    for item in keyL:
        print '\t'.join(['#'.join(i) for i in aDict[item]])
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


