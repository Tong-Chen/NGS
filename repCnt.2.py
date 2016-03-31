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
    if lensysargv != 5:
        print >>sys.stderr, "Find repeat in file"
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename.forC minlen \
maxlen[not include] repCnt' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    minlen = int(sys.argv[2])
    maxlen = int(sys.argv[3])
    repCnt = int(sys.argv[4])
    aDict = {}
    for line in open(sys.argv[1]):
        if line[0] != '>':
            line = line.strip()
            length = len(line)
            savedSet = set()
            stop = length - minlen + 1
            for start in range(stop):
                for width in range(minlen, maxlen):
                    end = start + width
                    if end > length:
                        break
                    seed = line[start:end]
                    if seed in savedSet:
                        continue
                    cntseed = line.count(seed)
                    if seed not in aDict:
                        aDict[seed] = cntseed
                    else:
                        aDict[seed] += cntseed
                    savedSet.add(seed)
                #------------------End one start--------------
            #------------End one sequence----------------------
        #----------------End one sequence------------------------
    #----------------End all sequences------------------------
    for key, value in aDict.items():
        if value >= repCnt:
            print "%s\t%s" % (key, value)
#----------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


