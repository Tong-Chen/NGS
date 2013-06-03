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
This is designed to bin transcriptome for coverage analysis.
Input file:
    chr18   36526168        36526187        NM_001081365_31.UTR5    0   +
    chr18   36552962        36553024        NM_001081365_31.UTR3    0   +
    chr17   26012444        26012474        NM_026686_32.UTR5       0   +
    chr17   26013942        26014108        NM_026686_32.UTR3       0   +
    chr18   36526187        36526395        NM_001081365_31.Coding_exon.1   0       +
    chr18   36552849        36552962        NM_001081365_31.Coding_exon.2   0       +
    chr17   26012474        26012671        NM_026686_32.Coding_exon.1  0       +
    chr17   26012904        26013067        NM_026686_32.Coding_exon.2  0       +
    chr17   26013164        26013224        NM_026686_32.Coding_exon.3  0       +
    chr17   26013439        26013507        NM_026686_32.Coding_exon.4  0       +
    chr17   26013605        26013684        NM_026686_32.Coding_exon.5  0       +
    chr17   26013894        26013942        NM_026686_32.Coding_exon.6  0       +

Output file:
    1.
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    desc = ""
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="Usually a bed file containing UTR5, \
Coding exon or UTR3.")
    parser.add_option("-l", "--length-of-bin", dest="len_bin",
        default=25, metavar=25, help="The length of bins you expected. \
Not exactly but roughly the numbner given here. The program will pick \
a suitable number to avoid very small length for last bin.")
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
    verbose = options.verbose
    debug = options.debug
    len_bin = int(options.len_bin)
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    for line in fh:
        lineL = line.split()
        start = int(lineL[1])
        end   = int(lineL[2])
        name  = lineL[3]
        len_r = end - start
        time  = len_r / len_bin
        #for regions with length smaller than expected bin
        if time == 0:
            lineL[3] = '__'.join([name, '1', '1'])
            print '\t'.join(lineL)
            continue
        real_len = len_r / time
        for i in range(time):
            lineL[1] = str(start + i * real_len)
            if i+1 == time:
                lineL[2] = str(end)
            else:
                lineL[2] = str(start + (i+1) * real_len)
            lineL[3] = '__'.join([name, str(time), str(i+1)])
            print '\t'.join(lineL)
        #----------The last one---------------
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



