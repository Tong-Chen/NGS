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
    if len(sys.argv) < 3:
        print >>sys.stderr, "The extended version of diffCom.py. You can \
    specify a colum to compare, all the colums will be output."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s singleColFile \
multipleColFile(tab separated) the special column(the column you want\
 to match, 1-based, default 1) the number of header lines[default 0] \
 auto_comp[complement missing value with 0, default no-complement. If \
 a parameter given, use it as the missed value]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------

    single   = sys.argv[1]
    multiple = sys.argv[2]
    if len(sys.argv) > 3:
        col = int(sys.argv[3]) - 1
    else:
        col = 0
    #---------------------------
    if len(sys.argv) > 4:
        header = int(sys.argv[4])
    else:
        header = 0
    #----------------------------
    if len(sys.argv) > 5:
        auto_comp = sys.argv[5]
    else:
        auto_comp = ''
    #-----------------------------------------
    singleL = []
    for line in open(single):
        singleL.append(line.strip())
    if auto_comp:
        singL_bak = singleL[:]
    #---------------------------------------------
    auto_linel = ''
    for line in open(multiple):
        if header:
            header -= 1
            print line,
            continue
        lineL = line.strip().split('\t')
        key = lineL[col]
        if key in singleL:
            print line
            if auto_comp:
                singL_bak.remove(key)
    #------missing-----------------------------
    if auto_comp and singL_bak:
        auto_linel = lineL[:]
        len_auto = len(auto_linel)
        for i in range(0, col):
            auto_linel[i] = auto_comp
        for i in range(col+1, len_auto):
            auto_linel[i] = auto_comp
        for key in singL_bak:
            auto_linel[col] = key
            print '\t'.join(auto_linel)
    #------------------------------------------
    #------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


