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

def readfile(file, start, end, size):
    aDict = {}
    aList = []
    for line in open(file):
        lineL = line.strip().split('\t')
        if lineL[5] == ',':
            continue
        valueL = [float(i) for i in lineL[5].split(',')[:-1]]
        tmpsum = sum(valueL[start:end])
        #duplicated lines are taken as one single
        aDict[line] = tmpsum
        aList.append(tmpsum)
    #---------------------------------------------
    bList = []
    aList.sort()
    diff = (aList[-1] - aList[0])/size
    for i in range(size):
        bList.append(aList[0] + (i+1) * diff)
    #make sure the biggest value included
    bList[size-1] = aList[-1] + 1
    #--------------------------------
    return aDict, bList
#--------------------------------------------------
def main():
    print >>sys.stderr, "This parse the CEAS dump file and separate \
regions into <given size> parts by the sum value of \
<specified region>."
    print >>sys.stderr, "Print the result to screen"
    lensysargv = len(sys.argv)
    if lensysargv != 5:
        print >>sys.stderr, 'Using python %s filename number of parts\
 start[1-based, included] end[1-based, included]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file = sys.argv[1]
    size = int(sys.argv[2])
    start = int(sys.argv[3])-1
    end = int(sys.argv[4])
    aDict, aList = readfile(file, start, end, size)
    fileList = [open(file+'.'+str(i+1),'w') for i in range(size)]
    for key, value in aDict.items():
        for i in range(size):
            if value < aList[i]:
                print >>fileList[i], key,
                break
        #----------finish one----------
    #-------fnish all-----close file------
    for fh in fileList:
        fh.close()
    #-------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


