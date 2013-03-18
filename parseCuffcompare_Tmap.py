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
        print >>sys.stderr, "This file is used to count number of \
various types Print the result to screen"
        print >>sys.stderr, 'Using python %s filename[- means \
sys.stdin]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    classCode = {
        'u':'Unknown_intergenic'
        'i':'Intron_derived'
        'j':'Novel_isoforms'
        'x':'Exon_antisense'
        'o':'Exonic_overlap'
        'c':'Contained'
        'p':'Run-on_frag'
        '=':'Known'
        'e':'Pre-mRNA'
        'r':'Repeats'
        's':'Intron_antisense'
        '.':'Multiple_classification'
        }


    file = sys.argv[1]
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #-------------------------
    header = 1
    for line in fh:
        if header:
            header -= 1
            continue
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


