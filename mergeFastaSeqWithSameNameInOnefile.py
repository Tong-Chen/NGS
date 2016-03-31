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

This is designed to merge fasta sequences with same names in one file.

Input:
>NM_001164147_26897@P_2625_1_2
GCCGTCACAGCTTCTTCGTTGGACCAGGGACCACCATATCTCTGCCCGGGGTGCATCAGGACGGTTCTGGA
>NM_001164147_26897@P_2625_1_2
CACAACCCAGCCGGGCACCTGTGA

Output:
>NM_001164147_26897@P_2625_1_2
GCCGTCACAGCTTCTTCGTTGGACCAGGGACCACCATATCTCTGCCCGGGGTGCATCAGGACGGTTCTGGACACAACCCAGCCGGGCACCTGTGA

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
        metavar="FILEIN", help="A fasta file, all strings all > will \
be treated as sequence names.")
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
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    for line in fh:
        if line[0] == '>':
            key = line[1:-1]
        else:
            if key not in aDict:
                aDict[key] = [line.strip()]
            else:
                aDict[key].append(line.strip())
    #-------------END reading file----------
    keyL = aDict.keys()
    keyL.sort()
    for key in keyL:
        print ">%s\n%s" % (key, ''.join(aDict[key]))
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



