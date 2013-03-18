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
    if lensysargv != 4:
        print >>sys.stderr, "This is used to bin regions for \
coverage computation. If <resolution> is positive value, \
<number> means the \
maximum bins for the longest region(s). If number of segments of one \
region less than this, overlapped regions will be used. This is used \
for regions \
with_same length or most with same length. \
If resolution is negative value, all regions will be split \
into <number> segments. if the segments less than <resolution> \
overlapped fragments will be produced."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s bedfile resolution \
number' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    resolution = int(sys.argv[2])
    number = int(sys.argv[3])
    mid = number / 2
    overlap = 1
    if resolution > 0:
        maxlen = resolution * number
    elif resolution < 0:
        minlen = resolution * number * (-1)
    #------------------------------------------
    for line in open(sys.argv[1]):
        lineL = line.strip().split("\t", 4)
        chr, name = lineL[0], lineL[3]
        other = ''
        if len(lineL) == 5:
            other = lineL[4]
        start = int(lineL[1])
        end   = int(lineL[2])
        tmpend = start
        mid = number / 2
        if resolution > 0:
            diff = maxlen + start - end
            diff = diff % number
            print "diff\t%d"  % diff
            if diff < number:
                shift = 1
            else:
                shift = diff / number + 1
            print "shift\t%d" % shift
            for i in range(number):
                if diff != 0 and (mid-diff/2.0)<=i and i<(mid+diff/2.0):
                    print "SHift--"
                    tmpstart = tmpend - shift
                    tmpend = tmpstart + resolution
                else:
                    tmpstart = tmpend
                    tmpend = tmpstart + resolution
                #-----------------------------------
                if other:
                    print "%s\t%d\t%d\t%s-%d\t%s" % \
                        (chr,tmpstart,tmpend,name,tmpstart,other)
                else:
                    print "%s\t%d\t%d\t%s-%d" % \
                        (chr,tmpstart,tmpend,name,tmpstart)
            #------------------------------------
        elif resolution < 0:
            minlen = resolution * number * (-1)
            diff = minlen - start + end
            for i in range(number):
               pass 
        #------------------------------------------------
    #----------------------------------------------------
#---------------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


