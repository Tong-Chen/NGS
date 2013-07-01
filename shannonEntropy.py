#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2013, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
'''
Functionla description

File format:
test_id TP22    TP23    TP16    eiNSC   sertoli
0610007C21Rik__XLOC_020042__=   7.14519 45.6536 46.4017 100.757 96.7285
0610007L01Rik__XLOC_020643__=   18.5009 7.2631  9.46834 11.1662 44.7087
0610007P08Rik__XLOC_006767__=   6.67136 5.67468 7.66136 4.43225 4.36321
0610007P14Rik__XLOC_006227__=   20.2475 26.611  40.441  58.4041 55.4284
0610007P22Rik__XLOC_011164__=   32.1198 55.6399 37.212  27.8336 18.7448
0610009B22Rik__XLOC_004475__=   8.67431 10.6717 13.0577 14.0596 22.4368
0610009D07Rik__XLOC_005329__=   46.6707 64.4448 70.8208 51.2847 43.1716
0610009O20Rik__XLOC_012558__=   8.14571 18.6023 17.8774 13.5515 13.2361
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
from math import log as ln

def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    desc = ""
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A tab separated file with the first \
column as names.")
    parser.add_option("-n", "--number_header", dest="header",
        metavar="1", default=1, help="A number to indicate the header \
lines in the file. Default 1 means only one header line.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    header = int(options.header)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    for line in fh:
        if header:
            header -= 1
            print "%s\tShannon_index" % line.rstrip()
            continue
        #-------------------------
        lineL = line.split()
        expr = [float(i) for i in lineL[1:]]
        if verbose:
            print >>sys.stderr, expr
        expr_sum = sum(expr)
        if verbose:
            print >>sys.stderr, expr_sum
        assert expr_sum != 0
        expr_R = [1.0 * i / expr_sum for i in expr]
        if verbose:
            print >>sys.stderr, expr_R
        expr_Log = []
        for i in expr_R:
            if i != 0:
                expr_Log.append(i*ln(i)/ln(2))
            else:
                expr_Log.append(i)
        if verbose:
            print >>sys.stderr, expr_Log
        shannon = -1 * sum(expr_Log)
        print "%s\t%s" % (line.strip(),str(shannon))
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    if verbose:
        print >>sys.stderr,\
            "--Successful %s" % strftime(timeformat, localtime())
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()



