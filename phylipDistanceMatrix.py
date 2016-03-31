#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
import sys
from math import sqrt
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"

def readDist(file, totalLine):
    dataL = []
    fh = open(file)
    species = int(fh.readline())
    assert totalLine % species == 0
    lineForSp = totalLine // species
    for j in range(species):
        tmpL = []
        for i in range(lineForSp):
            for item in fh.readline().split():
                tmpL.append(item)
        #-------------------------------------------------       
        dataL.append(tmpL[:])
    fh.close()
    #---------transfer to one way list-----------------------------------
    dataL1 = []
    for i in range(1, species):
        for j in range(1,i+1):
            dataL1.append(float(dataL[i][j]))
        #---------------------------------------
    return dataL1
#------------------------------------

def main():
    lensysargv = len(sys.argv)
    if lensysargv != 4:
        print >>sys.stderr, "This compute the correlation coefficient of \
    two things like proteins by supplied distance matrix[only lower \
    triangle is used]. The output of phylip protdist is the test matrix. \
    Usually it will ignore lines with only one element, such as the \
    sequence/species number."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s dist1 dist2 totalLine' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    totalLine = int(sys.argv[3]) - 1
    dataL1 = readDist(sys.argv[1], totalLine)
    #print sum(dataL1)
    dataL2 = readDist(sys.argv[2], totalLine)
    #print sum(dataL2)
    lenDL1 = len(dataL1)
    lenDL2 = len(dataL2)
    assert lenDL1 == lenDL2
    aevDL1 = sum(dataL1) / lenDL1
    aevDL2 = sum(dataL2) / lenDL2 
    #-------compute------------------------------------
    numerator = sum([(dataL1[i]-aevDL1) * (dataL2[i]-aevDL2) \
        for i in range(lenDL1)])
    denominator = sqrt(sum([(i-aevDL1) ** 2 for i in dataL1])) \
        * sqrt(sum([(i-aevDL2) ** 2 for i in dataL2])) 
    if denominator != 0:
        print numerator / denominator
#-------------------------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


