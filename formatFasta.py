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
    This is designed to transfer single-line FATSA or multiple lines
    FATSA into multiple-lines FASTA with each line no longer than
    given length.

Purpose:
    This is first designed for nafold since it demands the input file
    have FASTA sequence no more than 80 nt in each line.

Input Format:
    >test
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCC
    OR
    >test
    AAAAAAAAAAAAAA
    AAAAAAAAAAAAAA
    AAAAAAACCCCCCC
Ouutput format (when given length is 25):
    >test
    AAAAAAAAAAAAAAAAAAAAAAAAA
    AAAAAAAAAACCCCCCCCCCCCCCC
    CCCCC

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
        metavar="FILEIN", help="Any legal FASTA file.")
    parser.add_option("-l", "--max-len", dest="max_len",
        metavar="MAX_ALLOWED_LENGTH", help="An integer number \
to set the maximum allowed length for each sequence line.")
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
    max_len = int(options.max_len)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    seqL = []
    key = ''
    for line in fh:
        if line[0] == '>':
            if key and seqL:
                print key, 
                seq = ''.join(seqL)
                len_seq = len(seq)
                print '\n'.join(\
                    [seq[i:i+max_len] for i in range(0, len_seq, max_len)])
            #---------------------- 
            key = line
            seqL = []
        else:
            seqL.append(line.strip())
    #-------------END reading file----------
    if key and seqL:
        print key, 
        seq = ''.join(seqL)
        len_seq = len(seq)
        print '\n'.join(\
            [seq[i:i+max_len] for i in range(0, len_seq, max_len)])
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



