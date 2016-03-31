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
import re
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"

def output(seq, chr, lower):
    #bed file zero started, half closed, same with iterator.span()
    iterator = lower.finditer(seq)
    for match in iterator:
        print "%s\t%d\t%d" % (chr, match.start(), match.end())
#----------------------------------------
def main():
    print >>sys.stderr, "Get repeat sequences from genome downloaded \
from ucsc. lowercase letter represents repeats."
    print >>sys.stderr, "Print the result to stdout"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    lower = re.compile(r'[a-z]+')
    chr = ''
    adict = {}
    for line in open(sys.argv[1]):
        if line[0] == '>':
            if chr and chr != line[1:-1]:
                output(adict[chr], chr, lower)
                adict = {}
                chr = line[1:-1]
                adict[chr] = ''
            else:
                chr = line[1:-1]
                adict[chr] = ''
        else:
            adict[chr] += line.strip()
        #----------------------------------
    #-----------------------------------------
    output(adict[chr], chr, lower)
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    #startTime = localtime()
    #startTime = '-'.join([str(x) for x in startTime[:3]]) \
    #    + ' ' + ':'.join([str(x) for x in startTime[3:6]])
    main()
    endTime = strftime(timeformat, localtime())
    #endTime = localtime()
    #endTime = '-'.join([str(x) for x in endTime[:3]]) \
    #    + ' ' + ':'.join([str(x) for x in endTime[3:6]])
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


