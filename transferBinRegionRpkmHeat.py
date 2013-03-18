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
    if lensysargv != 2:
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    aDict = {}
    for line in open(sys.argv[1]):
        keyL, rpkm = line.strip().split("\t")
        mainKey, assistKey = keyL.split('-')
        if mainKey not in aDict:
            aDict[mainKey] = {}
        if assistKey not in aDict[mainKey]:
            aDict[mainKey][assistKey] = rpkm
        else:
            print >>sys.stderr, "duplicate"
        #-------------------------------------
    #-----------------------------------------
    mainKeyL = aDict.keys()
    mainKeyL.sort(key=lambda x:int(x))
    len_val = len(aDict[mainKeyL[0]].keys())
    print "Name\t%s" % ('\t'.join([str(i) for i in range(len_val)]))
    for mainKey in mainKeyL:
        assistKeyL = aDict[mainKey].keys()
        assistKeyL.sort(key=lambda x: int(x))
        print "%s\t%s" % (mainKey, \
            '\t'.join([aDict[mainKey][i] for i in assistKeyL]))
#-----------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


