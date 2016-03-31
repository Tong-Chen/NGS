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
    if lensysargv != 5:
        print >>sys.stderr, "Find repeat in file"
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename.forC minlen \
maxlen[not include] repCnt' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    minlen = int(sys.argv[2])
    maxlen = int(sys.argv[3])
    repCnt = int(sys.argv[4])
    fh = open(sys.argv[1])
    line = fh.readline()
    length = int(fh.readline())
    line = fh.readline()
    savedSet = set()
    stop = length - repCnt * minlen
    for start in range(stop):
        for length in range(minlen, maxlen):
            seed = line[start:start+length]
            if seed in savedSet:
                continue
            cntseed = line.count(seed)
            if cntseed >= repCnt:
                print "%s\t%d" % (seed, cntseed)
                savedSet.add(seed)
            else:
                break
        #------------------End one start--------------
    #------------End all----------------------
#----------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


