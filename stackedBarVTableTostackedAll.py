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
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    head = 1
    adict = {}
    for line in open(sys.argv[1]):
        if head:
            keyL = line.strip().split('\t')[1:]
            lenKeyL = len(keyL)
            for item in keyL:
                adict[item] = []
            head = 0
        else:
            lineL = line.strip().split('\t')
            sample = lineL[0]
            for i in range(lenKeyL):
                tmpkey = keyL[i]
                secondvalue = tmpkey.replace(' ','+')
                value = lineL[i+1]
                newitem = sample + '\t' + secondvalue + ':' + '\t' + \
                    value + '\t' + '=>\t' + value + '%' 
                adict[tmpkey].append(newitem) 
        #-----------------------------------------
    #------------------------------------------------
    i = 0
    for item in keyL:
        if i:
            print '=multi'
        print '\n'.join(adict[item])
        i = 1
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


