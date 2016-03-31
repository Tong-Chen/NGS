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

def combine(alpha, length):
    alist = []
    for item in alpha:
        
#-----------------end of combine------------------
def main():
    lensysargv = len(sys.argv)
    if lensysargv != 3:
        print >>sys.stderr, "Generate symbol combination from given \
alphabet and length."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s alphabet[AGCT] minlen[5]\
 maxlen[10]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    alpha = list(sys.argv[1])
    minlen = int(sys.argv[2])
    maxlen = int(sys.argv[3])
    for i in range(minlen, maxlen+1):
        combine(alpha, i)
#------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


