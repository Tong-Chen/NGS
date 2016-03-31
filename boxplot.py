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
Datafile
File 1
1
2
3
4
File 2
1
2
3
4

'''
import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
def main():
    lensysargv = len(sys.argv)
    if lensysargv < 8 or lensysargv % 2:
        print >>sys.stderr, "This uses boxplot.n in gplots package to do \
the boxplot of given files. Each file with one column. Filename and \
label must match."
        print >>sys.stderr, 'Using python %s boxplot.n[boxplot] \
prefix main xlab ylab filename1 2 ..., label1 2 ... ' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    func = sys.argv[1]
    prefix = sys.argv[2]
    main = sys.argv[3]
    xlab = sys.argv[4]
    ylab = sys.argv[5]
    mid = lensysargv/2+3
    fileL = sys.argv[6:mid]
    labelL = sys.argv[mid:]
    file = prefix+'.boxplot.r'
    #----------------------------------------
    fh = open(file, 'w')
    for i in range(mid-6):
        print >>fh, "%s <- read.table(file=\"%s\", header=FALSE)" % \
            (labelL[i], fileL[i])
    print >>fh, "data <- list(%s)" % ','.join([i+'='+i+"$V1" for i in labelL])
    print >>fh,"png(filename=\"%s\", width=960, height=960, res=150)" %\
            (prefix+'.boxplot.png')
    
    if func == "boxplot.n":
        print >>fh, '''
library(gplots)
boxplot.n(data, varwidth=TRUE, main=\"%s\", xlab=\"%s\", ylab=\"%s\")
''' % (main, xlab, ylab)
    elif func == 'boxplot':
        print >>fh, "boxplot(data, varwidth=TRUE, main=\"%s\",\
xlab=\"%s\", ylab=\"%s\")" % (main, xlab, ylab)
    #-----------------------------
    print >>fh, "dev.off()"
    fh.close()
    cmd = 'Rscript '+ file
    #print cmd
    os.system(cmd)
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


