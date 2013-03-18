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
#----------------format of each file------------------------
'''
small_bin_cnt (Only the first two columns are enough, the third one
unsed. **** Sort by the first column ***)
mm9_chr1_67108800       0       100
mm9_chr1_67108850       0       100
mm9_chr1_134217650      1       100
mm9_chr1_134217700      1       100
mm9_chr1_8388550        0       100
mm9_chr1_8388600        0       100
mm9_chr1_16777150       0       100

large_bin_cnt (Only the first two columns are enough, the third one unsed)
NR_038165_1     2
NR_038166_2     2
NM_027855_3     6
NM_001081394_4  3
NM_027854_5     3


'''

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
        print >>sys.stderr, 'Using python %s small_bin_cnt[- means \
sys.stdin,  ***Sort by first column***] large_bin_cnt \
small_bin_large_bin_relationship [***Sort by first column***]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file = sys.argv[1]
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    for line in fh:
        #here is your reading
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


