#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to generate a combination or permutation of a list of strings.

    For example: combination([1, 2, 3], 2) will generate '[(1,  2),  (1,  3),  (2,  3)]'
    permutations([1, 2, 3], 2]) will generate '[(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]'
'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from itertools import permutations, combinations, combinations_with_replacement

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0


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
        metavar="FILEIN", help="A file with specified column containing strings to be combined.")
    parser.add_option("-H", "--header", dest="header",
        type="int", default=0, help="Number of header lines to be skipped. Default no header line.")
    parser.add_option("-c", "--col-number", dest="col",
        type="int", default=1, help="Default 1 meaning the first column.")
    parser.add_option("-o", "--operation", dest="operation",
        type="choice", choices=["permutations", "combinations", "combinations_with_replacement"], 
        help="The operation to be done, including <combinations_with_replacement>, <combinations>, <permutations>")
    parser.add_option("-n", "--items-number-to-comb", dest="num",
        type="int", help="Number of items used for combination. Optional only for <permutations>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    header = options.header
    col  = options.col - 1
    operation = options.operation
    num  = options.num
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    itemL = [line.strip().split('\t')[col] for line in fh.readlines()]
    itemL = itemL[header:]
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    opD = {'permutations': permutations, 'combinations': combinations, 
           'combinations_with_replacement': combinations_with_replacement}
    if num:
        resultL = opD[operation](itemL, num)
    else:
        if operation == 'permutations':
            resultL = permutations(itemL)
        else:
            print >>sys.stderr, "Please specify <--items-number-to-comb> with a non-0 interger."
            return 1
    #-----output-----------
    for result in resultL:
        print '\t'.join(result)
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
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
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


