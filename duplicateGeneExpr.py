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
    This is designed to deal with expression matrix with duplicate gene names.
'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
import numpy as np

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
        metavar="FILEIN", help="A gene expression matrix.")
    parser.add_option("-m", "--method", dest="method",
        help="For genes with duplicate expression values, choose the suitable method to keep only one expression value for each gene. Currently, <sum (expression summary)>, <mean (mean expression), <var (gene with most variable expressions)>, <random (randomly keep one)> are supported.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def dealwithExprLL(exprLL, method):
    '''
    exprLL = [[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]]
    '''
    if method == 'mean':
        return [str(i) for i in np.mean(exprLL, 0)]
    elif method == 'sum':
        return [str(i) for i in np.sum(exprLL, 0)]
    elif method == 'var':
        varL = np.var(exprLL, 1)
        posL = np.argsort(varL)
        return [str(i) for i in exprLL[posL[-1]]]
#------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    method  = options.method
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
            print line,
            continue
        #---------------------------
        lineL = line.strip().split()
        gene  = lineL[0]
        exprL = [float(i) for i in lineL[1:]]
        if gene not in aDict:
            aDict[gene] = []
        aDict[gene].append(exprL)
    #---------------------------------------
    for gene, exprLL in aDict.items():
        if len(exprLL) == 1 or method == 'random':
            exprL = [str(i) for i in exprLL[0]]
        else:
            exprL = dealwithExprLL(exprLL, method)
        print "{}\t{}".format(gene,'\t'.join(exprL))
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


