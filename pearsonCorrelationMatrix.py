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

    Output would be <filein>.<method>.xls.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from scipy.stats.stats import pearsonr
import pandas as pd
from numpy import log2

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
        metavar="FILEIN", help="A matrix with first row as header and \
first column as row names. The results will be output to files.")
    parser.add_option("-r", "--row", dest="row_correlation",
        default=False, action="store_true", help="Default perform column correlation. \Specify to compute row correlation.")
    parser.add_option("-l", "--log2-transform", dest="log2_trans",
        default=False, action="store_true", 
        help="Perform log2 transform before computing pearson correlation value. A pseudo count 1 will be added before log2 transform.")
    parser.add_option("-m", "--method", dest="method",
        default="pearson", help="pearson (default), kendall, spearman")
    #parser.add_option("-s", "--scale", dest="scale",
    #    default=False, action="store_true" help="Specify to scale data before correlation analysis. For row_correlation, scale will be perfomed on columns. For column correlation,  scale will be performed on rows.")
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
    row_correlation = options.row_correlation
    log2_trans = options.log2_trans
    method = options.method
    output = file + '.'+method+'.xls'
    #fh_out = open(output, 'w')
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    matrix = pd.read_table(file, header=0, index_col=0)
    if log2_trans:
        matrix = log2(matrix+1)
    if row_correlation:
        matrix = matrix.T
    corr = matrix.corr(method=method)
    corr.index.name = method
    corr.to_csv(output, sep=b'\t')
    #-------------END reading file----------
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


