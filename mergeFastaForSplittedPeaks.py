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

This is designed to paste sequences together.

InputFile:
    1.fa
    >NM_001081956_3160__1__S__2
    TGGAAGAAGAGAGGAAAAAAGAAAAGAAGAGAGAAGAGGAAGTAGCTGAAAG
    >NM_001081956_3160__1__S__1
    AAAAAGAAAAGATGAAGTAAGGAAAGCCCAAGAAAAGAGAAAGAAGGCCAG
    >NM_010127_21849__4__S__1
    GTGCAGCCTATCCAGCCGACACAAGCCGTGCCCCAGCCTGCAGT
    >NM_010127_21849__4__S__2
    TGCCCAGGGCCAGGTGATCGCAACCCTAGCCA
    2.bed (only  the forth and sixth columns are used.)
	chr15	100411268	100411382	NM_010127_21849__4__S__1	2.6	-
	chr15	100413658	100413758	NM_010127_21849__4__S__2	2.6	-
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
    parser.add_option("-i", "--fasta", dest="fasta",
        metavar="fasta", help="A fasta file with all sequences in one \
line is needed here.")
    parser.add_option("-b", "--bed", dest="bed",
        metavar="bed", help="The bed file used to get the fasta \
sequence.")
    parser.add_option("-s", "--sep", dest="sep",
        metavar="a symbol", default="__", help="The separtor used to separate \
FASTA sequence names and bed line names to get common part and\
 diff part. Sequences will be merged together if they have \
 same common part in their names. The diff part is a number,  which \
 will be used to get the right order of those sequences. \
 Default '__' for our example sample.")
    parser.add_option("-c", "--com_index", dest="com_index",
        default="1,2,3", help="The part of string you want as common \
part (1-based). Default '1,2,3' for our example.")
    parser.add_option("-d", "--diff_index", dest="diff_index",
        default="4", help="The part of string you want as common \
part (1-based). Default '4' for our example.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    (options, args) = parser.parse_args(argv[1:])
    assert options.fasta != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    fasta = options.fasta
    bed   = options.bed
    sep   = options.sep
    com_index = [int(i)-1 for i in options.com_index.split(',')]
    diff_index = int(options.diff_index)-1
    verbose = options.verbose
    #-----------------------------------
    fastaD = {}
    if fasta == '-':
        fh = sys.stdin
    else:
        fh = open(fasta)
    #--------------------------------
    for line in fh:
        if line[0] == '>':
            key = line[1:-1]
        else:
            fastaD[key] = line.strip()
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    bedD = {}
    #---------Read bed---------------
    for line in open(bed):
        lineL = line.strip().split()
        name = lineL[3]
        strand = lineL[5]
        if strand == '+':
            strand = 1
        elif strand == '-':
            strand = -1
        else:
            print >>sys.stderr, "Wrong strand %s" % strand
            sys.exit(1)
        nameL = name.split(sep)
        com_name = sep.join([nameL[i] for i in com_index])
        diff_name = int(nameL[diff_index]) * strand
        if com_name not in bedD:
            bedD[com_name] = {}
        assert diff_name not in bedD[com_name]
        bedD[com_name][diff_name] = name
    #-------------------------------------------------------
    #----------sort--------------------
    for key, valueD in bedD.items():
        diffKeyL = valueD.keys()
        diffKeyL.sort()
        print '>%s' % key
        print ''.join(\
            [fastaD[valueD[diffkey]] for diffkey in diffKeyL])
    #----------output--------------------------

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



