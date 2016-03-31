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
        metavar="FILEIN", help="miRNA list, only the ones before first \
blank need to be suppled here.")
    parser.add_option("-s", "--mature-seq", dest="seq",
        metavar="FILEIN", help="Fatsa file get from miRBase. The name \
before the first blank will be used as the name for fasta seq. \
The '>' is excluded.")
    parser.add_option("-r", "--reverse_complement", dest="rc",
        metavar="0/1", default=0, help="If 1, the target sequence of miRNA \
seed will be returned instead of seed. Default 0.")
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
    seq  = options.seq
    rc   = options.rc
    verbose = options.verbose
    debug = options.debug
    dict = {'A':'U','G':'C','U':'A','C':'G','a':'u','g':'c','u':'a','c':'g'}
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    miRD = {}
    for line in fh:
        miRD[line.strip()] = 1
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    for line in open(seq):
        if line[0] == '>':
            output = 0
            key = line[1:].split()[0]
            if key in miRD:
                output = 1
                miRD.pop(key)
        elif output == 1:
            seed = line[1:8]
            if rc:
                target = ''.join([dict[i] for i in seed[::-1]])
                print "%s\t%s" % (key, target)
            else:
                print "%s\t%s" % (key, seed)
    #------------------------------------------
    if miRD:
        print >>sys.stderr, "Unmatched mir %s" % '\t'.join(miRD.keys())
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



