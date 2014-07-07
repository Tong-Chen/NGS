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
    This is designed to transfer predicted miRNA target genome
    coordinates in bed format.
    The input file cound be all miRanda output or just the lines begin
    with only one '>' since only these lines will be used.
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
        metavar="FILEIN", help="The input file cound be all \
miRanda output or just the lines begin with only one '>' \
since only these lines will be used.")
    parser.add_option("-b", "--bed", dest="bed",
         help="The bed file contains the coordinates of \
sequences used for target prediction. Six column file with \
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
        if line[0] == '>' and line[1] != '>':
            lineL = line[1:].split()
            key = lineL[1]
            name = ':'.join([key, lineL[0]])
            if name in nameDict:
                nameDict[name] += 1
            else:
                nameDict[name] = 1
            tmpL = bedDict[key]
            strand = tmpL[5]
            if strand == '+':
                start = int(tmpL[1]) + int(lineL[6]) - 1
                end   = int(tmpL[1]) + int(lineL[7])
            elif strand == '-':
                start = int(tmpL[2]) - int(lineL[7])
                end   = int(tmpL[2]) - int(lineL[6]) + 1
            if nameDict[name] == 1:
                print "%s\t%s\t%s\t%s\t%s\t%s" % (tmpL[0], str(start),
                        str(end), name, lineL[2], strand)
            else:
                print "%s\t%s\t%s\t%s:%d\t%s\t%s" % (tmpL[0], str(start),
                        str(end), name, nameDict[name], lineL[2], strand)

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



