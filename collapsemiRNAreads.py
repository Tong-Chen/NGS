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
    This script is designed to transfer smRNA-Seq reads in
    'reads<tab>count' format to fasta format for 
    <quantifier.pl> in mirDeep2.

    Use <mapper.pl> from mirDeep2 to transferring FATSA 
    or FASTQ format to <quantifier.pl> input,
    like <mapper.pl reads_normal_fasta -c -m -s reads_collapsed_fa>
    or   <mapper.pl reads_normal_fastq -e -m -s reads_collapsed_fa>


    Sample file for 'read<tab>count' format:
    #ID_REF = 
    #VALUE = number of times sequenced
    ID_REF  VALUE
    ATAGAATATAACCTTTGCGTGT  6
    GTTTTATTCGGCAAAGAG      10


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
        metavar="FILEIN", help="Normal read-count format file. Usually \
        can be directly downloaded from NCBI.")
    parser.add_option("-H", "--header", dest="header",
        help="The number of header lines in given file. \
        No default value.")
    parser.add_option("-s", "--species", dest="species",
        default='seq', help="A three alphabet letter to indicate the \
        name of the species. Optional default 'seq' to be consistent \
        with default value in <mapper.pl>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    assert options.header != None, "Please specify number of header lines"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    header = int(options.header)
    prefix = options.species
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    num = 1
    for line in fh:
        if header:
            header -= 1
            continue
        #------------------
        seq, count = line.split()
        print ">%s_%s_x%s\n%s" % (prefix,num,count, seq)
        num += 1
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



