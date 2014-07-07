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
    Order FASTA files alphabetically or according to a list given in a
    file.
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
        metavar="FILEIN", help="A FASTA file.")
    parser.add_option("-o", "--order-file", dest="order_file",
        metavar="ORDER_FILE", help="Only given when you want to sort \
FASTA file as this order. This file contains the names of FATSA \
sequence in files given to -i. ")
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
    order = options.order_file
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    orderL = []
    if order:
        orderL = [line.strip() for line in open(order)]
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    add = 0
    for line in fh:
        if line[0] == '>':
            add = 0
            key = line[1:].strip()
            if (not orderL) or (orderL and key in orderL):
                aDict[key] = []
                add = 1
            #----------------------------------------
        else:
            if add:
                aDict[key].append(line)
    #-------------END reading file----------
    if orderL:
        for key in orderL:
            print ">%s" % key
            print ''.join(aDict[key]),
    else:
        keyL = aDict.keys()
        keyL.sort()
        for key in keyL:
            itemL = aDict[key]
            print ">%s\n%s" % (key, ''.join(itemL).rstrip('\n'))
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



