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

    This is designed to extract miranda hits from miranda output. If
    only for this purpose, you can just use [grep -B 16 '^>>' output].
    This program accepts a file contains the miRNA-gene pair to label
    those verified target pairs and extract them.
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
        metavar="FILEIN", help="The output of miranda")
    parser.add_option("-f", "--pair-file", dest="pair",
        metavar="PAIR-FILE", help="A two column file contains miRNA \
and target gene pairs with miRNAs at the first column.")
    parser.add_option("-o", "--output-unverified-pair", dest="out_all",
        metavar="OUT-ALL", default=0, help="Default FALSE meaning only \
output verified pair. Accept any string represents TRUE to output all \
predicted pairs after a cutting line. Also when this is TRUE, a \
parameter to -f may not be needed.")
    parser.add_option("-s", "--sta", dest="sta",
        default=0, help="Count the miRNA-target gene pairs which are \
both verified and predicted.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readMiranda(fh):
    aDict = {}
    while 1:
        line = fh.readline()
        while not line.startswith("Performing Scan:"):
            line = fh.readline()
            if line.startswith("Scan Complete"):
                return aDict
        lineL = line.split()
        key1 = lineL[2]
        key2 = lineL[4]
        if key1 not in aDict:
            aDict[key1] = {}
        if key2 not in aDict[key1]:
            aDict[key1][key2] = [line]
        #---------------------------------
        line = fh.readline()
        while not line.startswith("Complete"):
            aDict[key1][key2].append(line)
            line = fh.readline()
        #----END one pair------------
    #--------END all pairs-----------
#----------------------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    miranda = options.filein
    pair    = options.pair
    out_all = options.out_all
    sta     = options.sta
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if miranda == '-':
        fh = sys.stdin
    else:
        fh = open(miranda)
    #--------------------------------
    aDict = readMiranda(fh)
    if out_all == 0 or pair:
        for line in open(pair):
            mir, gene = line.split()
            tmp = ''.join(aDict[mir][gene])
            if tmp.find("No Hits Found above Threshold") != -1:
                print 'v-',''.join(aDict[mir][gene])
            else:
                print 'v+',''.join(aDict[mir][gene])
            if out_all:
                aDict[mir].pop(gene)
        #-------------END reading file----------
    if out_all:
        print '--------------------------'
        for mir, itemD in aDict.items():
            for gene, itemL in itemD.items():
                tmp = ''.join(itemL)
                if tmp.find("No Hits Found above Threshold") == -1:
                    print 'uv+',''.join(aDict[mir][gene])
                
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



