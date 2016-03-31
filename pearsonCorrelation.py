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
    This uses scipy.stats.stats.pearsonr to compute pearson
    correlation coefficient.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from scipy.stats.stats import pearsonr

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
        metavar="FILEIN", help="A matrix. The correlation \
will be computed for rows. The results will be output to files.")
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
    output = file + '.scipy.pearson'
    fh_out = open(output, 'w')
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    header = 1
    aDict = {}
    for line in fh:
        if header:
            header -= 1
            continue
        #--------------------------
        lineL = line.strip().split('\t')
        key = lineL[0]
        valueL = [float(i) for i in lineL[1:]]
        aDict[key] = valueL
    #-----END reading-----------------------
    print >>fh_out, "Item1\tItem2\tPearson_correlation\tp"
    keyL = aDict.keys()
    len_keyL = len(keyL)
    for i in range(len_keyL-1):
        key1 = keyL[i]
        value1 = aDict[key1]
        for j in range(i+1, len_keyL):
            key2 = keyL[j]
            value2 = aDict[key2]
            cor, p = pearsonr(value1, value2)
            print >>fh_out, "\t".join([key1, key2, str(cor), str(p)])
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    fh_out.close()
    cmd = ["multipleTest.sh -f", output, "-p holm"]
    os.system(' '.join(cmd))
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


