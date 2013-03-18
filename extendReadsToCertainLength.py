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

def read_chromosome(file):
    return dict([line.strip().split('\t') for line in open(file)])

#-------------------------------------------
def main():
    print >>sys.stderr, "Extend reads to given length."
    print >>sys.stderr, "Print the result to screen"
    lensysargv = len(sys.argv)
    if lensysargv != 4:
        print >>sys.stderr, 'Using python %s bedfile \
genomeFile(chromosome size) finalLength' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    aDict = read_chromosome(sys.argv[2])
    size = int(sys.argv[3])
    for line in open(sys.argv[1]):
        lineL = line.strip().split('\t')
        if lineL[5] == '-':
            lineL[1] = int(lineL[2]) - size
            if lineL[1] < 0:
                lineL[1] = 0
            lineL[1] = str(lineL[1])
        elif lineL[5] == '+':
            lineL[2] = int(lineL[1]) + size
            if lineL[2] > aDict[lineL[0]]:
                lineL[2] = aDict[lineL[0]]
            lineL[2] = str(lineL[2])
        #-----------------------------------------
        print '\t'.join(lineL)
    #---------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


