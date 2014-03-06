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
    This is designed to get random lines from a file.
'''

import sys
import os
from random import randint
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
        metavar="FILEIN", help="")
    parser.add_option("-n", "--number-of-output-lines", dest="num",
        metavar="An Integer", help="A number to indicate the number \
of lines you want to extract from the file. Header lines were \
not counted and files with less than these lines will give \
a warning.")
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
    num = int(options.num)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
        print >>sys.stderr, "No STDIN allowed."
        sys.exit(1)
    else:
        fh = open(file)
    #--------------------------------
    if verbose:
        print >>sys.stderr, "Total # of wanted lines is ", num
    header = 1
    count = 0
    for line in fh:
        if header:
            header -= 1
            continue
        #--------------------
        count += 1
    #-------------END reading file----------
    if verbose:
        print >>sys.stderr, "Total number of lines is ", count
    if file != '-':
        fh.close()
    if count <= num:
        for line in open(file):
            print line,
        if count < num:
            print >>sys.stderr, "Warning: no enough line for output."
    else:
        random = count - num #the maximum allowed depleted lines
        if verbose:
            print >>sys.stderr, \
                "The maximum allowed depleted lines is", random
        real_output = 0
        header = 1
        for line in open(file):
            if header:
                header -= 1
                print line,
                continue
            #--------------------------
            if random != 0:
                select = randint(0, 1)
                if verbose:
                    print >>sys.stderr, \
                        "The fate of current line is", select
                if select:
                    print line,
                    real_output += 1
                else:
                    random -= 1
            else:
                    print line,
                    real_output += 1
            #---------------------------
            if real_output == num:
                break
    #----close file handle for files-----
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



