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
    print >>sys.stderr, "Count specific character[s] by line."
    print >>sys.stderr, "Print the result to screen."
    lensysargv = len(sys.argv)
    if lensysargv < 5:
        print >>sys.stderr, 'Using python %s filename[colums \
separated by tab] \
nameColumn[1-based] countColumn[1-based] countChar[s] Num[\
accept an integer, specially 1 means count the terms]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file=sys.argv[1]
    nameCol=int(sys.argv[2]) - 1
    countCol = int(sys.argv[3]) - 1
    countChar = sys.argv[4]
    num = int(sys.argv[5]) if lensysargv>5 else 0
    #-----------------------------------------
    for line in open(file):
        lineL = line.strip().split("\t")
        print "%s\t%d" % (lineL[nameCol], \
            lineL[countCol].count(countChar)+num)
    #----------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


