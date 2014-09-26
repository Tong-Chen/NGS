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
    This is designed to select differentially expressed genes when no
    repeats are available.
    Normally, it accepts a matrix containing gene expression in
    multiple samples and compares gene expression differences between 
    each two samples. 
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
        metavar="FILEIN", help="An expression matrix file with each \
line represents gene expression and each column represents \
special sample. No blank is allowed in each column. A header line \
with sample names are required.")
    parser.add_option("-f", "--fold", dest="fc",
        default=2, help="The fold-change used to specify DE genes. \
A float or int number is required, with <2> as default.")
    parser.add_option("-l", "--log2", dest="log2",
        default=1, help="This parameter is used to indicate if \
expression value is in log2 format. Default <TRUE> means the \
program assumes the expression data is in log2 format and will \
transfer it to normal representation. If a FALSE value given, \
no additional preprocess needed.")
    parser.add_option("-s", "--simple", dest="simple",
        default=0, help="When a value (like 1) representing TRUE \
is given, only compare the numbers no matter fold-change. \
Higher expressed genes and lower expressed genes will be both output.")
    parser.add_option("-j", "--jitter", dest="jitter",
        default=1, help="Add a value to all expression values \
to avoid divide to 0. Default 1.")
    parser.add_option("-p", "--prefix", dest="prefix",
        help="The prefix for output files. Default filename \
given to -i.")
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
    fc   = float(options.fc)
    fc_r = 1 / fc
    simple = options.simple
    log2 = options.log2
    jitter = float(options.jitter)
    prefix = options.prefix
    if prefix == None:
        prefix = file
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    head = 1
    fhDict = {}
    for line in fh:
        if head:
            sampL = line.split()[1:]
            lenSampL = len(sampL)
            for i in range(lenSampL-1):
                if i not in fhDict:
                    fhDict[i] = {}
                for j in range(i+1, lenSampL):
                    fhDict[i][j] = {}
                    fileU = '.'.join([prefix, sampL[i], 'up', sampL[j]]) 
                    fileD = '.'.join([prefix, sampL[i], 'dw', sampL[j]]) 
                    fhDict[i][j]['up'] = open(fileU, 'w')
                    fhDict[i][j]['dw'] = open(fileD, 'w')
            #--------------------------------------
            head -= 1
            continue
        #------------------------------------
        lineL = line.split()
        gene = lineL[0]
        if log2 and (not simple):
            exprL = [2**float(i) + jitter for i in lineL[1:]]
        else:
            exprL = [float(i) + jitter for i in lineL[1:]]
        #---------------------------------
        for i in range(lenSampL-1):
            expr_i = exprL[i]
            for j in range(i+1, lenSampL):
                expr_j = exprL[j]
                time = expr_i / expr_j
                if not simple:
                    if time >= fc:
                        print >>fhDict[i][j]['up'], gene
                    elif time <= fc_r:
                        print >>fhDict[i][j]['up'], gene
                else:
                    if time > 1:
                        print >>fhDict[i][j]['up'], gene
                    elif time < 1:
                        print >>fhDict[i][j]['up'], gene
            #----------------------------------------------
        #--------------------------------------------------
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    for key, valueD in fhDict.items():
        for secKey, secValueD in valueD.items():
            for thirdV in secValueD.values():
                thirdV.close()
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



