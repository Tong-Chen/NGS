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
    if lensysargv < 2:
        print >>sys.stderr, "This program is used for sorting one file \
    according to another file."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s file_need_to_sort(-)  \
file_need_to_sort_header_lines file_ordered col_ordered' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file_unsorted = sys.argv[1]
    header = int(sys.argv[2])
    file_ordered = sys.argv[3]
    col_ordered  = int(sys.argv[4]) - 1
    #------------------------------------
    sortL = [line.split()[col_ordered] for line in open(file_ordered)]

    headL = []
    aDict = {}
    if file_unsorted == '-':
        fh = sys.stdin
    else:
        fh = open(file_unsorted)
    #-------------------------
    for line in fh:
        if header:
            headL.append(line)
            header -= 1
        #---------------------
        key = line.split()[0]
        if key not in aDict:
            aDict[key] = line
        else:
            print >>sys.stderr, "Duplicate key %s" % key
    #------END for -------------
    if file_unsorted != '-':
        fh.close()
    #-------close --------------
    if headL:
        print ''.join(headL),
    for key in sortL:
        if key in aDict:
            print aDict[key],
        else:
            print >>sys.stderr, "Unexisted key %s" % key

#----------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


