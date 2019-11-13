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
    This is designed to compute the sample score for given gene set.

    Normally,  log2 transformed normalized gene expression data for given gene set will be averaged. Then these score will be centered and divided by their standard deviation.
    If multiple sets of genes are supplied, then these scores for each cell will be centered and divided by their standard deviation.

Ref:
    We then averaged the normalized expression levels (log2(TPM+1)) of the genes in each gene-set to define the phase-specific scores of each cell. These scores were then subjected to two normalization steps. First,  for each phase,  the scores were centered and divided by their standard deviation. Second,  the normalized scores of each cell were centered and normalized.

EXPR_mat:
    ID  samp1   samp2   samp3   samp4
    A   1000    2000    1500    1700
    B   1000    2000    1500    1700
    C   1000    2000    1500    1700
    D   1000    2000    1500    1700

Gene_set (one set or multiple sets)
The first column is set names, remaining columns are IDs (which should match the 
first column of EXPR_mat). Tab separated
    Set1    A   B   C    
    Set2    C   D
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
import pandas as pd
from numpy import log2

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

def fprint(content):
    """ 
    This is a Google style docs.

    Args:
        param1(str): this is the first param
        param2(int, optional): this is a second param
            
    Returns:
        bool: This is a description of what is returned
            
    Raises:
        KeyError: raises an exception))
    """
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
        metavar="FILEIN", help="EXPR mat")
    parser.add_option("-o", "--output-file", dest="output",
        help="Name for output file")
    parser.add_option("-g", "--gene-set", dest="gene_set",
        help="A file containing gene sets with format specified above.")
    parser.add_option("-I", "--ignore-missing-gene", dest="ignore_missing_gene",
        action="store_true", default=False,  help="ignore_missing_gene. Default False.")
    parser.add_option("-l", "--log2-transform", dest="log2_trans",
        default=0, type='float', 
        help="Whether or not perform log2 transform for expression value. \
For expression values which have already transformed, ignore this parameter. \
Otherwise a non-zero value should be given here, which would be added to all \
expr values and then log2-transform performed on each expr value.")
    parser.add_option("-r", "--raw", dest="raw",
        action="store_true", default=False,  help="Output raw score. Default False.")
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
    output = options.output
    gene_set = options.gene_set
    gene_setL = [line.strip().split('\t', 1) for line in open(gene_set)]
    log2_trans = options.log2_trans
    raw = options.raw
    ignore_missing_gene = options.ignore_missing_gene
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if debug:
        print >>sys.stderr, gene_setD
    expr_mat = pd.read_table(file, header=0, index_col=0)
    expr_mat.index = expr_mat.index.map(str)
    if log2_trans:
        expr_mat = log2(expr_mat+log2_trans)
    scoreL = []
    for key, geneL in gene_setL:
        geneL = geneL.split('\t')
        tmp_expr_mat = expr_mat[expr_mat.index.isin(geneL)]
        geneL = set(geneL)
        miss_gene = geneL.difference(set(tmp_expr_mat.index))
        if miss_gene:
            print >>sys.stderr, \
                "Following genes in {} are missing: {}".format(key, miss_gene)
            if not ignore_missing_gene:
                sys.exit(1)
        #------------------------------------
                
        tmp_expr_mat_mean = tmp_expr_mat.mean()
        tmp_expr_mat_mean = tmp_expr_mat_mean.to_frame(name=key)
        if debug:
            tmp_expr_mat_mean.to_csv(key+'.tmp', sep=b'\t')
        if not raw:
            tmp_expr_mat_mean = (tmp_expr_mat_mean-tmp_expr_mat_mean.mean())/tmp_expr_mat_mean.std()
        scoreL.append(tmp_expr_mat_mean)
    #------------------------------------
    resultM = pd.concat(scoreL,  axis=1, join="outer")
    if not raw and len(gene_setL)>1:
       resultM = resultM.T
       mean_tmp = resultM.mean()
       std_tmp  = resultM.std()
       resultM_norm = (resultM-mean_tmp) / std_tmp
       resultM = resultM_norm.T
    resultM.index.name = 'ID'
    resultM.to_csv(output, sep=b"\t")

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


