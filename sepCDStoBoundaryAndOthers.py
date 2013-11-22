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

This program is written to get the boundary part and other parts for
all blocks. The first boundary and last boundary is not considered.


Input file format:

chr7	52823164	52823749	NR_038165_1.Coding_exon.1	0	-
chr7	52829782	52829892	NR_038165_1.Coding_exon.2	0	-
chr7	52829977	52830147	NR_038165_1.Coding_exon.3	0	-
chr7	52830496	52830546	NR_038165_1.Coding_exon.4	0	-
chr7	52823164	52823749	NR_038166_2.Coding_exon.1	0	-
chr7	52826355	52826562	NR_038166_2.Coding_exon.2	0	-
chr7	52829782	52829892	NR_038166_2.Coding_exon.3	0	-
chr7	52829977	52830147	NR_038166_2.Coding_exon.4	0	-
chr7	52830496	52830546	NR_038166_2.Coding_exon.5	0	-
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
        metavar="FILEIN", help="A bed file like given sample")
    parser.add_option("-s", "--size-for-boundary", dest="size",
        metavar="30", default=30, help="A number to represent the \
length of the boundary. All inner exons with length larger than 60 \
will get two boundaries. The first and last exon will get one inner \
boundary.")
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
    size = int(options.size)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    for line in fh:
        
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



