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
    This is designed to transfer relative postions to absolute
    positions.

    Assume we have a file contains many sequences and we searched some
    short strings and got their relative positions. We also have a bed
    file contains the coordinates of all sequences. Then we can use
    this program to get the coordinates of all appeared short strings.

    RelativePositionFile (The second and third column represent the start
    and end position of short strings of the forth column in sequences
    nameed as the first column. These positions are 0-started,
    start_positions instead of end_positions are included. The forth
    column is not needed):

    P_29109 108     113     GGACT
    P_24697 47      52      GGACT
    P_3675_2        1       6       GGACA
    P_3675_2        17      22      GGACA
    P_3675_2        115     120     GAACA
    P_3675_2        127     132     GGACA
    P_24691 16      21      GGACC
    
    SequenceCoordinateFile:

    chr8    87385225        87385300        P_29019 1.47323  +
    chr11   57872225        57872400        P_3675_2    3.2184  +
    chr6    125314601       125314675       P_24697 2.31  -
    chr6    125235985       125236060       P_24691 2.6  -

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
    usages = "%prog -i file -b bedfile"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="The input file should contain \
at least three columns with the first column as the name in \
the forth column given to -b.")
    parser.add_option("-b", "--bed", dest="bed",
         help="The bed file contains the coordinates of \
sequences used for short motif search. Six column file with \
strand information needed.")
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
    bed  = options.bed
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    bedDict = {}
    for line in open(bed):
        lineL = line.split()
        key = lineL[3]
        assert key not in bedDict, "Duplicate keys %s" % key
        bedDict[key] = lineL
    #------------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    nameDict = {}
    for line in fh:
        lineL = line.split()
        len_lineL = len(lineL)
        name = lineL[0]
        if name in nameDict:
            nameDict[name] += 1
        else:
            nameDict[name] = 1
        tmpL = bedDict[name]
        strand = tmpL[5]
        if strand == '+':
            start = int(tmpL[1]) + int(lineL[1])
            end   = int(tmpL[1]) + int(lineL[2])
        elif strand == '-':
            start = int(tmpL[2]) - int(lineL[2])
            end   = int(tmpL[2]) - int(lineL[1])
        if nameDict[name] == 1:
            if len_lineL < 4:
                print "%s\t%s\t%s\t%s\t%s\t%s" % (tmpL[0], str(start),
                        str(end), name, tmpL[4], strand)
            else:
                print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (tmpL[0], str(start),
                        str(end), name, tmpL[4], strand,
                        '\t'.join(lineL[3:]))
        else:
            if len_lineL < 4:
                print "%s\t%s\t%s\t%s:%d\t%s\t%s" % (tmpL[0], str(start),
                    str(end), name, nameDict[name], tmpL[4], strand)
            else:
                print "%s\t%s\t%s\t%s:%d\t%s\t%s\t%s" % \
                    (tmpL[0], str(start),
                    str(end), name, nameDict[name], tmpL[4], strand, 
                    '\t'.join(lineL[3:]))
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



