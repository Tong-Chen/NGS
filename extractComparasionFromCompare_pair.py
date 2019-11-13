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
    This is designed to extract a subset of comparison from complete
    combination of comparison. 

Input_file 1 (full comparison, given to -i)

TR42261|c0_g1   T0_vs_T2.T0-UP
TR74130|c2_g1   T0_vs_T4.T0-UP
TR17401|c0_g1   T0_vs_T12.T0-UP
TR65645|c1_g2   T0_vs_T24.T0-UP
TR18148|c0_g1   T0_vs_T48.T0-UP
TR84404|c0_g1   T0_vs_T2.T2-UP
TR74130|c2_g1   T0_vs_T4.T4-UP
TR17401|c0_g1   T0_vs_T12.T12-UP
TR65645|c1_g2   T0_vs_T24.T24-UP
TR18148|c0_g1   T0_vs_T48.T48-UP

Input_file 2 (compare_pair, given to -c)

T0      T12
T0      T24
T0      T48


'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from bs4 import BeautifulSoup

#reload(sys)
#sys.setdefaultencoding('utf8')

#from multiprocessing.dummy import Pool as ThreadPool

debug = 0

def fprint(content):
    print json_dumps(content,indent=1)

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
        metavar="FILEIN", help="Two columns file as described above.")
    parser.add_option("-c", "--compare-pair", dest="compare_pair",
        metavar="COMPARE-PAIR", help="Two columns file as described above.")
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
    compare_pair = options.compare_pair
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    pairD = {}
    for line in open(compare_pair):
        comp1, comp2 = line.split()
        comp = '_'.join([comp1, 'vs', comp2])
        comp1 = comp + '.' + comp1 + '-UP'
        comp2 = comp + '.' + comp2 + '-UP'
        pairD[comp1] = 1
        pairD[comp2] = 1
    #------------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    head = 0
    for line in fh:
        if head:
            print line,
            head -= 1
            continue
        #-----------------
        lineL = line.rstrip().rsplit('\t', 1)
        key = lineL[1]
        if key in pairD:
            print line,
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
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


