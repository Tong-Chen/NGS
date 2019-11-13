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
    This is designed to label samples by their expression values.


Input data:

gene	col44A7_90	col44A7_96	col44A2_6	col44A8_21	col44A1_50	col44A8_23	col44A8_25
LGR5    1   1   2   2   3   3   2
SOX2    1   1   2   2   3   3   2
PROM1   1   1   2   2   3   3   2
POU5F1   1   1   2   2   3   3   2
CD44    1   1   2   2   3   3   2

Program:
    labelSamplesByGeneExpression.py -i input.file -g 'LGR5;PROM1;CD44'

Output data:

    <table format>

    samp	LGR5	PROM1	CD44
    col44A7_90	neg	neg	neg
    col44A8_21	pos	pos	pos
    col44A1_50	pos	pos	pos
    col44A8_23	pos	pos	pos
    col44A8_25	pos	pos	pos
    col44A7_96	neg	neg	neg
    col44A2_6	pos	pos	pos

    <two_col format>

    col44A7_90	LGR5.neg_PROM1.neg_CD44.neg
    col44A8_21	LGR5.pos_PROM1.pos_CD44.pos
    col44A1_50	LGR5.pos_PROM1.pos_CD44.pos
    col44A8_23	LGR5.pos_PROM1.pos_CD44.pos
    col44A8_25	LGR5.pos_PROM1.pos_CD44.pos
    col44A7_96	LGR5.neg_PROM1.neg_CD44.neg
    col44A2_6	LGR5.pos_PROM1.pos_CD44.pos

'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

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
        metavar="FILEIN", help="An expression matrix")
    parser.add_option("-g", "--gene", dest="gene",
        help="';' separated multiple genes used for labeling samples.")
    parser.add_option("-c", "--classifier", dest="classifier",
        default="average", 
        help="Accept <average> (default) or <median> to split samples into two groups, higher than average value and lower than average value. Accept <quantile> to split samples into foru groups (lower than first quantile, first quantile to median, median to third quantile, larger than third quantile)")
    parser.add_option("-F", "--format", dest="format", default="two_col", 
        help="Specify output format, <table> format containing full information for each sample. <two_col> format including two columns with classification information integrated. See above for examples.")
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
    file  = options.filein
    gene  = options.gene
    geneL = gene.split(';')
    classifier = options.classifier
    format = options.format
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    geneL2 = geneL[:]
    header = 1
    for line in fh:
        lineL = line.split()
        if header:
            aDict = {}
            headerL = lineL[1:]
            for i in headerL:
                aDict[i] = {}
            header -= 1
            continue
        #---------------------------
        #print >>sys.stderr, geneL2
        #print >>sys.stderr, geneL
        gene = lineL[0]
        if gene in geneL:
            exprL = [float(i) for i in lineL[1:]]
            ave = sum(exprL)/len(exprL)
            for samp, expr in zip(headerL, exprL):
                if classifier == "average":
                    if expr >= ave:
                        aDict[samp][gene] = 'pos'
                    else:
                        aDict[samp][gene] = 'neg'
                #------------------------------
            #------------------------------
            geneL.remove(gene)
            if not geneL:
                break
        #--------------------------------------
    #-------------END reading file----------

    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    #print >>sys.stderr, geneL2
    if format == "table":
        print "samp\t{}".format('\t'.join(geneL2))
        for samp, sampD in aDict.items():
            tmpL = [sampD[gene] for gene in geneL2]
            tmpL.insert(0, samp)
            print "\t".join(tmpL)
    elif format == "two_col":
        print "samp\tconditions"
        for samp, sampD in aDict.items():
            tmpL = [gene+'.'+sampD[gene] for gene in geneL2]
            tmp = '_'.join(tmpL)
            print "{}\t{}".format(samp, tmp)
        
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


