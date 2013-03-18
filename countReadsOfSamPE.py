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
from os import system
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
def main():
    print >>sys.stderr, "This is used to count the mapped reads \
in SAM/BAM output of paired-end reads."
    print >>sys.stderr, "Print the result to file"
    lensysargv = len(sys.argv)
    if lensysargv != 3:
        print >>sys.stderr, 'Using python %s filename suffix[filetype, \
SAM or BAM]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file = sys.argv[1]
    output = file + '.sta'
    com =  file + " | cut -f 2 | sort | uniq -c \
| sed 's/^  *//' | awk 'BEGIN{OFS=\"\\t\" ;FS=\" \"}{print \
$2,$1}' > " + output

    type = sys.argv[2].upper()
    if type == "BAM":
        cmd="samtools view -X " + com
    elif type == "SAM":
        cmd= "cat " + com
    else:
        print >>sys.stderr, "Wrong type, need to be BAM or SAM"
    #----------------------------------------------------------------
    for line in open(output):
        if line.find('p'):
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


