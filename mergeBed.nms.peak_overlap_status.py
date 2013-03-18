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
'''
Input file

chr1    3027303 3028464 MACS_1422_1;MACS_513_1;MACS_17-3-15_1
chr1    3111487 3112435 MACS_17-3-15_2;MACS_1422_2;MACS_F1_6-1_1;MACS_513_2;MACS_233_1
chr1    3137130 3138441 MACS_R1_1;MACS_17-3-15_3
chr1    3434754 3435420 MACS_F1_MEF_1
chr1    3480242 3480872 MACS_F1_6-1_2
chr1    3504603 3505910 MACS_F1_MEF_2

For this input, targetColumn is 4. separtor_in_targetColumn is ';'.
separtor_for_identity is '_'. identity is the second element.
'''
def main():
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, "This is used to summary peak overlap \
status among multiple samples. It will output the number of peaks \
coocupied by 1, 2, 3 , ... samples. Default parameters are suitable for \
mergeBed output using -nms parameter only.Print the result to screen"
        print >>sys.stderr, 'Using python %s filename \
targetColumn[1-based, default 4] separtor_in_targetColumn[ default ;] \
separtor_for_identity[default _] position_of_identity[default 2, \
1-based]' % sys.argv[0]
        print >>sys.stderr, "====Sample line-------"
        print >>sys.stderr, "chr1    3137130 3138441    MACS_R1_1;MACS_17-3-15_3"
        sys.exit(0)
    #-----------------------------------
    file = sys.argv[1]
    if lensysargv > 2:
        targetCol = int(sys.argv[2]) - 1
    else:
        targetCol = 3
    if lensysargv > 3:
        sepForCol = sys.argv[3]
    else:
        sepForCol = ';'
    if lensysargv > 4:
        sepForIden = sys.argv[4]
    else:
        sepForIden = '_'
    if lensysargv > 5:
        posOfIden = int(sys.argv[5]) - 1
    else:
        posOfIden = 1
    #----------------------------------------------
    aDict = {}
    for line in open(file):
        tmpL = set([i.split(sepForIden)[posOfIden] for i in \
                (line.split("\t")[targetCol]).split(sepForCol)])
        num = len(tmpL)
        if num not in aDict:
            aDict[num] = 1
        else:
            aDict[num] += 1
        #-----------------------------------------------------
    #-----------------------------------------------------
    keyL = aDict.keys()
    keyL.sort()
    for key in keyL:
        print "%d\t%d" % (key, aDict[key])
#--------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


