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
        metavar="FILEIN", help="")
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
    firstD = {}
    secondD = {}
    for line in fh:
        first,second = line.split()
        if first not in firstD:
            firstD[first]= set()
        firstD[first].add(second)
        if second not in secondD:
            secondD[second] = set()
        secondD[second].add(first)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    keyL = firstD.keys()
    rmSet = set()
    for key in keyL:
        if key in firstD:
            valueS = firstD[key]
            len1 = len(valueS)
            if len1 < 1:
                continue
            rmkey = ''
            while rmkey != key and valueS:
                for key70 in valueS:
                    if key70 in firstD:
                        tmpLen = len(firstD[key70])
                        if tmpLen > len1:
                            rmkey = key70
                #-----------------------------
                if rmkey == '':
                    rmkey = key
                else:
                    valueS.discard(rmkey)
                rmSet.add(rmkey)
                firstD.pop(rmkey)
                for key85 in firstD.keys():
                    firstD[key85].discard(rmkey)
        #------------------------------------
    #-------------------------------------------
    print '\n'.join(rmSet)
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



