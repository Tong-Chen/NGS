#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division  #, with_statement
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

def readfile(file, key, aDict):
    assert key not in aDict
    aDict[key] = {}
    for line in open(file):
        lineL = line.strip().split('\t')
        if len(lineL) == 1:
            keyI = lineL[0]
            vI = '1'
        else:
            keyI = lineL[0]
            vI = lineL[1]
        assert keyI not in aDict[key]
        aDict[key][keyI] = vI
    #------------------------------
#----------------------------------

def main():
    lensysargv = len(sys.argv)
    if lensysargv < 5 or not lensysargv % 2: 
        #lensysargv must be odd number
        print >>sys.stderr, "This program was first designed to \
collapse TE peak results of different samples into one file to \
compare the differences. If the item in main file did not exist in one \
subfile, 0 would be supplied as the value. If only contains one \
column in subfile, 1 will be given as value for iterms exist in \
mainfile and subfile."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s mainfile columnOfMainfile\
[1-based] \
 subfile1,...[the first column for match column in main, the second \
column is output, 1-based. If only one column, 1 is the output value \
default no header lines] subfilesymbol1,... ' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    mainfile=sys.argv[1]
    maincol = int(sys.argv[2]) - 1
    #mid = (lensysargv-3) / 2 + 3
    mid = (lensysargv+3) / 2
    subfiles=sys.argv[3:mid]
    subfilesymbol=sys.argv[mid:]
    keyL = [line.strip().split('\t')[maincol] for line in open(mainfile)]
    keyL = set(keyL)
    keyL = list(keyL)
    keyL.sort()
    aDict = {}
    lensubfiles = len(subfiles)
    for i in range(lensubfiles):
        #readfile(file, key, aDict):
        readfile(subfiles[i], subfilesymbol[i], aDict)
    #--------------------------------------------------
    print "Sample\t%s" % '\t'.join(subfilesymbol)
    for keyTe in keyL:
        tmpL = [keyTe]
        for i in range(lensubfiles):
            tmpDict = aDict[subfilesymbol[i]] 
            if keyTe in tmpDict:
                tmpL.append(tmpDict[keyTe])
            else:
                tmpL.append('0')
        #--------------------------------------
        print "\t".join(tmpL)
    #-------------------------------------
#----------------end main------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


