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
if False:
    print "This program does not work under python 3, \
run in python 2.x."

import sys
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
def main():
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename[- means \
sys.stdin] forParseGTF(given anything to turn on this option, optional)' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file = sys.argv[1]
    if lensysargv == 3:
        keepRows = ['exon', 'start_codon', 'stop_codon'] 
    else:
        keepRows = []

    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    aDict = {}
    for line in fh:
        lineL = line.split("\t")
        if (keepRows and lineL[2] in keepRows) or ((not keepRows) and lineL[2] != 'CDS'):
            firstKey = lineL[0] + '.'.join(lineL[8].split("; ")[0:2])
            if firstKey not in aDict:
                aDict[firstKey] = {}
            secondKey = (int(lineL[3]), int(lineL[4]), lineL[2])
            assert secondKey not in aDict[firstKey], \
                (firstKey, secondKey)
            aDict[firstKey][secondKey] = line
        #here is your reading
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    firstKeyL = aDict.keys()
    firstKeyL.sort()
    for firstKey in firstKeyL:
        secondKeyL = aDict[firstKey].keys()
        secondKeyL.sort()
        for secondKey in secondKeyL:
            print aDict[firstKey][secondKey],
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


