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

def find_all(astr, sub):
    start = 0
    while True:
        start = astr.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)
#-------------------------------------------------------



def main():
    print >>sys.stderr, "Print the result to screen"
    lensysargv = len(sys.argv)
    if lensysargv < 1:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    str = "CGACGCCCCTCC"
    print str
    print list(find_all(str, 'C'))
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


