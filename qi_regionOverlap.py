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
desc = '''
Functional description:
    This is designed to compute the overlap between two regions.

chr10   4892012 4892031 A15     CTCF    MA0139.1        chr10   4892012 4892038
chr17   4316664 4316683 A25     CTCF    MA0139.1        chr17   4316645 4316685

Only column 2,3,8,9 will be used, other columns will be outputted as
original.

These regions are default in Bed format,  [start, end) start is
included while end is excluded. 

'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A file contains at least 4 columns \
indicating the start and end of two regions.")
    parser.add_option("-r", "--region1", dest="region1",
        help="The column number of start and end for the first \
region. Like '2,3' in example for the first region.")
    parser.add_option("-R", "--region2", dest="region2",
        help="The column number of start and end for the second \
region. Like '8,9' in example for the second region.")
    parser.add_option("-m", "--minimum-allowed-overlap", dest="mo",
        help="Output only lines with overlapped regions no \
less than given number here.")
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
    r1 = [int(i.strip()) for i in options.region1.split(',')]
    r2 = [int(i.strip()) for i in options.region2.split(',')]
    mi_over = int(options.mo)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    for line in fh:
        lineL = line.split()
        r1_start = int(lineL[r1[0]-1])
        r1_end   = int(lineL[r1[1]-1])
        r2_start = int(lineL[r2[0]-1])
        r2_end   = int(lineL[r2[1]-1])
        r1_set = set()
        for i in range(r1_start,  r1_end):
            r1_set.add(i)
        r2_set = set()
        #-------------------------------------
        for i in range(r2_start,  r2_end):
            r2_set.add(i)
        #---------------------------------------------
        inter = len(r1_set.intersection(r2_set))
        if inter >= mi_over:
            print line,
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



