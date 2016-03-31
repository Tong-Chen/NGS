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
'''
Region  1422    2737
mm9_chr10_10000000      2       0
mm9_chr10_100000000     7       20
mm9_chr10_100001000     8       16
mm9_chr10_100003000     1       1
mm9_chr10_100005000     33      19
mm9_chr10_100007000     11      9
mm9_chr10_100008000     14      20
mm9_chr10_100009000     19      29
mm9_chr10_10001000      6       3
mm9_chr10_100010000     8       12
mm9_chr10_100011000     14      7
'''

import sys
from fisher import pvalue
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
def main():
    lensysargv = len(sys.argv)
    if lensysargv < 4:
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename total_first_col \
total_second_col two_tail[left_tail, right_tail] head[number of lines\
 needs to skip, default 1] divide_a_value[used when your count larger \
than 4294967296]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file = sys.argv[1]
    total_first_col = int(sys.argv[2])
    total_second_col = int(sys.argv[3])
    #print total_first_col,  total_second_col
    if lensysargv > 4:
        tail = sys.argv[4]
    else:
        tail = "two_tail"
    if lensysargv > 5:
        head = int(sys.argv[5])
    else:
        head = 1
    if lensysargv > 6:
        scale = int(sys.argv[6])
        total_first_col = total_first_col/scale
        total_second_col = total_second_col/scale
    else:
        scale = 1
    #-----------------------------------------------
    for line in open(file):
        line = line.rstrip()
        if head:
            print "%s\t%s" % (line, 'p')
            head -= 1
            continue
        #----------------------------------
        lineL = line.split()
        q = int(lineL[1]) / scale
        m = int(lineL[2]) / scale
        if q == 0 and m == 0:
            continue
        p = pvalue(q, m, total_first_col-q, total_second_col-m)
        if tail == 'two_tail':
            print "%s\t%s" % (line, p.two_tail)
        elif tail == 'left_tail':
            print "%s\t%s" % (line, p.left_tail)
        elif tail == 'right_tail':
            print "%s\t%s" % (line, p.right_tail)
#-----------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


