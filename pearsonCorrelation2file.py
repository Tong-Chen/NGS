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
from math import log
from statsmodels.stats.multitest import multipletests

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
        metavar="FILEIN", help="A matrix file. The correlation \
will be computed for each row in this file to all rows in file given to <-j>.")
    parser.add_option("-j", "--input-file2", dest="filein2",
        metavar="FILEIN", help="A matrix file. The correlation \
will be computed for each row in this file to all rows in file given to <-i>.")
    parser.add_option("-L", "--label", dest="label",
        metavar="FILEIN", help="Label for files given to <-i> and <-j>. \
Coma separated two strings, order matters. Like <file_label1,file_label2>. \
Default <Item1,Item2>.")
    parser.add_option("-l", "--log2-transform", dest="log2_trans",
        default=0, 
        help="Default 0. Accept a non-zero number to initiate log2 transform and this value will be added to all values to avoid log2 transform of 0.")
    parser.add_option("-c", "--compare-pair", dest="pairfile",
        help="Only compute correlation for given pairs. The first column relative to file given to <-i>. Second column relative to file given to <-j>. Other columsn will be ignored.")
    parser.add_option("-t", "--threshold", dest="threshold",
        default=0.6, type="float", 
        help="R square larger than 0.6 (default).")
    parser.add_option("-p", "--positive-cor-only", dest="pos_cor",
        default=False, action="store_true", 
        help="Specify to get only positive correlations.")
    parser.add_option("-o", "--output-prefix", dest="output",
        help="Prefix of output files")
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
    file2 = options.filein2
    pairfile = options.pairfile
    log2_trans = options.log2_trans
    threshold = options.threshold
    pos_cor = options.pos_cor

    if log2_trans:
        log2 = log(2)
        log2_trans = float(log2_trans)
    output = options.output + '.scipy.pearson.xls'
    if options.label:
        labelL = [i.strip() for i in options.label.split(',')]
    else:
        labelL = ["Item1", "Item2"]
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
        if log2_trans:
            valueL = [log(float(i)+log2_trans)/log2 for i in lineL[1:]]
        else:
            valueL = [float(i) for i in lineL[1:]]
        aDict[key] = valueL
    #-----END reading-----------------------
    bDict = {}
    header = 1
    for line in open(file2):
        if header:
            header -= 1
            continue
        #--------------------------
        lineL = line.strip().split('\t')
        key = lineL[0]
        if log2_trans:
            valueL = [log(float(i)+log2_trans)/log2 for i in lineL[1:]]
        else:
            valueL = [float(i) for i in lineL[1:]]
        bDict[key] = valueL
    #-----END reading-----------------------
    pairD = {}
    for line in open(pairfile):
        first, second = line.strip().split('\t')[:2]
        if first in aDict and second in bDict:
            pairD[(first, second)] = 1
    #------------------------------------------
    #--------------------------    
    labelL.extend(["Pearson_correlation", "R_square", "pvalue", "fdr"])

    print >>fh_out, "\t".join(labelL)
    keyL = aDict.keys()
    keyL.sort()
    key2L = bDict.keys()
    key2L.sort()
    resultL = []
    pL = []
    for key1 in keyL:
        value1 = aDict[key1]
        for key2 in key2L:
            if (key1, key2) not in pairD:
                continue
            value2 = bDict[key2]
            if debug:
                print >>sys.stderr, key1, value1
                print >>sys.stderr, key2, value2
            cor, p = pearsonr(value1, value2)
            if pos_cor and cor < 0:
                continue
            r2 = cor ** 2
            if r2 >= threshold:
                resultL.append([key1, key2, str(cor), str(r2), str(p), 1])
                pL.append(p)
            #print >>fh_out, "\t".join([key1, key2, str(cor), str(p)])
    #-------------END reading file----------
    p_adjL = multipletests(pL, method="fdr_bh")[1]
    for eachtmpL, p_adj in zip(resultL, p_adjL):
        eachtmpL[-1] = format(p_adj, '0.2E')
    resultL.sort(key=lambda x: float(x[4]))
    #----close file handle for files-----
    print >>fh_out, '\n'.join(['\t'.join(tmpL) for tmpL in resultL])
    if file != '-':
        fh.close()
    fh_out.close()
    #cmd = ["multipleTest.sh -f", output, "-p holm"]
    #os.system(' '.join(cmd))
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


