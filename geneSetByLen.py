#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division #, with_statement
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

def output(file, aDict, keyL, rangeL, grpNum):
    fhL = [open(file+'.'+str(i+1),'w') for i in range(grpNum)]
    for value in keyL:
        for i in range(grpNum):
            if i == 4:
                print >>fhL[i], '\n'.join(aDict[value]) 
                break
            elif value <= rangeL[i]:
                print >>fhL[i], '\n'.join(aDict[value]) 
                break
    #-----------------------------------------------------------
    [fh.close() for fh in fhL]
#----------------------------------------


def main():
    print >>sys.stderr, "This separate genes into given number sets by \
its length by two methods. One way, get the minimum and maximum number, \
use the quotient as the ranges for each gene set (diff). The other way is \
get unique gene length in a list, then separate the list to the given \
number set with each set has even number of different lengths (evenL).\
The second way is seperating gene sets into even number sets, \
size of sets is roughly quotient of totoal gene number divide group \
number. If genes with same length are separated into diffeent sets, \
they will be moved to the second set instead of first (evenN). \
This \
also can deal with expression data or any other data with the \
first column as label, the second column as value."
    print >>sys.stderr, "Print the result to screen"
    lensysargv = len(sys.argv)
    if lensysargv != 4:
        print >>sys.stderr, 'Using python %s filename[two column file] \
grpNum diff[evenL, evenN]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file = sys.argv[1]
    grpNum = int(sys.argv[2])
    type = sys.argv[3]
    aDict = {}
    count = 0 #This is only used for type 'evenN'
    for line in open(file):
        count += 1 #This is only used for type 'evenN'
        value, key = line.strip().split()
        key = float(key)
        if key in aDict:
            aDict[key].append(value)
        else:
            aDict[key] = [value]
    #--------------------------------------
    keyL = aDict.keys()
    keyL.sort()
    rangeL=[]
    #print aDict
    #print keyL
    #------diff--------------------
    if type == 'diff':
        minimum = keyL[0]
        maximum = keyL[-1]
        size = (maximum-minimum) / grpNum
        for i in range(grpNum-1):
            rangeL.append(keyL[0] + (i+1)*size)
    #------evenL--------------------
    elif type == 'evenL':
        length = int(len(keyL) / grpNum)
        for i in range(grpNum-1):
            rangeL.append(keyL[(i+1)*length-1])
    #------evenN--------------------
    elif type == 'evenN':
        length = int(count / grpNum)
        tmpRangeL = [(i+1)*length for i in range(grpNum-1)]
        #print tmpRangeL
        i = 0
        tmpCount = 0
        for value in keyL:
            if i == 4: break
            tmpCount += len(aDict[value])
            if tmpCount >= tmpRangeL[i]:
                tmpCount = tmpRangeL[i]
                rangeL.append(value)
                i = i + 1
        #-----------------------------------
        #print rangeL
    #--------------------------------------------
    output(file, aDict, keyL, rangeL, grpNum)
#------------------------------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


