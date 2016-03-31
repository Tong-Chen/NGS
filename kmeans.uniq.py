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
    print >>sys.stderr, "Make the specificed column uniq by add .i. \
This is written for R which read.table needs unique rowname. \
Print the result to screen"
    lensys = len(sys.argv)
    if lensys < 2 :
        print >>sys.stderr, 'Using python %s filename header[default 1] \
rowname[1-based column number, default 1]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    header = int(sys.argv[2])   if lensys == 3 else 1
    row    = int(sys.argv[3])-1 if lensys == 4 else 0
    aDict = {}
    for line in open(sys.argv[1]):
        if header:
            header -= 1
            continue
        #------------------------------------------
        lineL = line.strip().split()
        key = lineL[row]
        if key in aDict:
            aDict[key].append(line)
        else:
            aDict[key] = [line]
    #-----------------------------------------------
    for key, valueL in aDict.items():
        j = 1
        for line in valueL:
            newKey = key + '.' + str(j)
            j += 1
            line.replace(key, newKey)
            print line, 
    #----------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


