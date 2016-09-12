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
    This is designed for two purposes:

    1. Transfer FATSA files with sequence 
       in multiple line format to single line fasta format.
    2. Transfer any type of FASTA files to table format with 
       the first column containing sequence names and the
       second column containing sequences.

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
        metavar="FILEIN", help="Multiple line FATSA.")
    parser.add_option("-o", "--outformat", dest="ofmt",
        default='single_line_fasta', 
        help="Default `single_line_fasta` to transfer to format 1, \
accept `table` to transfer to format 2.")
    parser.add_option("-s", "--nameSep", dest="nameSep",
        default="haha_no_thisSeP_cT", help="Default using originally \
full name as new name. Accept a separtor to extract part of original \
name as new name.")
    parser.add_option("-p", "--position-index", dest="pos_index",
        default=1, help="Default <1> representing get the first part of original name \
as new name. Accept other numbers to extract corresponding parts.")
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
    ofmt = options.ofmt
    nameSep = options.nameSep
    pos_index = int(options.pos_index)-1
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    tmpL = []
    key = ''
    for line in fh:
        if line[0] == '>':
            if key:
                seq = ''.join(tmpL)
                if ofmt == "single_line_fasta":
                    print "%s\n%s" % (key, seq)
                elif ofmt == "table":
                    print "%s\t%s" % (key[1:], seq)
                tmpL = []
            key = line.strip().split(nameSep)[pos_index]
        else:
            tmpL.append(line.strip())
    #-------------END reading file----------
    if key:
        seq = ''.join(tmpL)
        if ofmt == "single_line_fasta":
            print "%s\n%s" % (key, seq)
        elif ofmt == "table":
            print "%s\t%s" % (key[1:], seq)
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



