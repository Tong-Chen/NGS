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
    This is designed to get the expression of genes from
    freeCompare.py output.

Input file:

    1. freeCompare.py or freeCompareWithName.py output
    CG11790 SC_11bA.SS._vs_.SC_dw__comm__SC_can.SS._vs_.SC_dw
    CG34439 SC_11bA.SS._vs_.SC_dw__comm__SC_can.SS._vs_.SC_dw
    Cpr66D  SC_11bA.SS._vs_.SC_dw__comm__SC_can.SS._vs_.SC_dw
    CG34437 SC_11bA.SS._vs_.SC_dw__comm__SC_can.SS._vs_.SC_dw
    CG31404 SC_11bA.SS._vs_.SC_dw__comm__SC_can.SS._vs_.SC_dw
    CG43235 SC_11bA.SS._vs_.SC_dw__comm__SC_can.SS._vs_.SC_dw
    CG31406 SC_11bA.SS._vs_.SC_dw__comm__SC_can.SS._vs_.SC_dw
    CG34433 SC_11bA.SS._vs_.SC_dw__comm__SC_can.SS._vs_.SC_dw

    2. Expr matrix output from DESeq2.sh
    
    3. Annotation file (with the first column matching the first
    column of file 1 and 2 as listed above)

    4. sampleFile
        Samp    conditions
        SC_1    SC
        SC_2    SC
        SC_3    SC
        SC_11bA.SS_1    SC_11bA.SS
        SC_11bA.SS_2    SC_11bA.SS
        SC_can.SS_1     SC_can.SS
        SC_can.SS_2     SC_can.SS
        SG_1    SG
        SG_2    SG
        SG_3    SG
        SG_bam_1        SG_bam
        SG_bam_2        SG_bam


'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
debug = 0
from math import log as ln

def fprint(content):
    print json_dumps(content,indent=1)

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = '''%prog -i free.all -e expr_matrix -s sampleFile -a anno

%prog -i free.all -e expr_matrix -s sampleFile -a anno -I all.desp
'''
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="The output of <freeCompare.py> or \
<freeCompareWithName.py>. Normally <*.free.all>. \
For the output of <freeCompareWithName.py>, \
an additional parameter to -I may be needed")
    parser.add_option("-I", "--input-file-desp", dest="filein_desp",
        metavar="FILEIN-DESP", help="Another output of  \
<freeCompareWithName.py> in format like <*.free.all.desp>.")
    parser.add_option("-e", "--expr-matrix", dest="expr",
        metavar="EXPR-MATRIX", help="The output of DESeq2.sh.")
    parser.add_option("-s", "--sampleFile", dest="samp",
        metavar="sampleFile", help="The sampleFile given to DESeq2.sh.")
    parser.add_option("-r", "--remove-string", dest="remove",
        metavar="REMOVE-STRING", default="union,comm,special", 
        help="Default <union, comm, special>, unused currently. \
All strings that are not sample in expr matrix skipped automatically.")
    parser.add_option("-R", "--replace-string", dest="replace",
        metavar="REPLACE-STRING", default="_up,_dw", help="Default <_up, _dw>")
    parser.add_option("-S", "--separtor", dest="sep",
        metavar="SEPARTOR", default="__,._vs_.", help="Default <__, ._vs_.>")
    parser.add_option("-a", "--anno", dest="anno",
        help="Annotation file")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def readSamp(samp):
    header = 1
    sampD = {}
    for line in open(samp):
        if header:
            header -= 1
            continue
        #-----------------------
        value, key = line.strip().split()
        if key not in sampD:
            sampD[key] = [value]
        else:
            sampD[key].append(value)
    #-------------------------------
    return sampD
#----------------------------------

def parseKey(key, sep, remove, rep, sampD):
    if debug:
        print >>sys.stderr, "Original", key
    for i in rep:
        key = key.replace(i, '')
    #---------replace--------------
    if debug:
        print >>sys.stderr, "Replace", rep
        print >>sys.stderr, "After replace", key
    keyL = [key]
    for i in sep:
        tmpKeyL = []
        for key in keyL:
            tmpKeyL.extend(key.split(i))
        keyL = tmpKeyL[:]
    #----------sep-------------------
    if debug:
        print >>sys.stderr, "Separtor", sep
        print >>sys.stderr, "After sep", keyL
    tmpKeyL = []
    for key in keyL:
        sampL = sampD.get(key, "unknown")
        if sampL != "unknown":
            tmpKeyL.extend(sampL)
    #------------remove--------------------
    if debug:
        print >>sys.stderr, "sampD", sampD
        print >>sys.stderr, "Final return", tmpKeyL
    return tmpKeyL
#------------------------------------------
def readMatrix(expr):
    if debug:
        start_count = 1
    header = 1
    matrixD = {}
    for line in open(expr):
        lineL = line.strip().split('\t')
        key = lineL[0]
        if header:
            header -= 1
            key = 'head'
            headerL = lineL
            matrixD[key] = headerL[1:]
            continue
        #-----------------------
        assert key not in matrixD, "Duplicate %s" % key
        matrixD[key] = {}
        lenLineL = len(lineL)
        for i in range(1, lenLineL):
            matrixD[key][headerL[i]] = lineL[i]
        if debug:
            if start_count < 10:
                print >>sys.stderr, "Key\tvalue", key, lineL
                #print >>sys.stderr, matrixD
            start_count += 1
    #---------------------------------
    return matrixD
#------------------------------------------
def readAnno(anno):
    annoD = {}
    if not anno:
        return annoD
    header = 1
    for line in open(anno):
        lineL = line.split('\t', 1)
        key = lineL[0]
        if header:
            key = 'head'
            header -= 1
        assert key not in annoD, key
        annoD[key] = line.strip()
    return annoD
#-------------------------------------

#def output(aDict, matrixD, prefix, annoD):
#    '''
#    aDict = {'compare': {
#                'samp':[]
#                'id':[]}   
#            }
#    matrixD = {'head':[],
#             'gene':{'samp':expr, 'samp2':expr2}}
#    '''
#    for key, valueD in aDict.items():
#        output = prefix+key+".results"
#        if annoD:
#            output_anno = prefix+key+".anno.xls"
#            anno_fh = open(output_anno, 'w')
#        fh = open(output, 'w')
#        sampleL = valueD['samp']
#        idL = valueD['id']
#        headerL = matrixD['head']
#        existL = [samp for samp in headerL if samp in sampleL]
#        print >>fh, "%s\t%s" % ('gene', '\t'.join(existL))
#        #print >>sys.stderr, matrixD['CG11790']
#        if annoD:
#            print >>anno_fh, "%s\t%s\t%s" % \
#                ('gene', '\t'.join(existL), annoD['head'])
#        for id in idL:
#            print >>fh, "%s\t%s" % (id, '\t'.join(\
#                [matrixD[id][samp] for samp in existL]))
#            if annoD:
#                print >>anno_fh, "%s\t%s\t%s" % (id, '\t'.join(\
#                    [matrixD[id][samp] for samp in existL]),
#                    annoD.get(id, ""))
#            #------END one id of each file------------
#        #----------END each item-----------------------
#        fh.close()
#        anno_fh.close()
##---------------------------------------------

def computeShannon(aList, plus=1):
    #print aList
    if len(aList) < 2:
        print >>sys.stderr, "You may need to specify \
-I parameter if the program stops."
    expr = [float(i)+1 for i in aList]
    expr_sum = sum(expr)
    assert expr_sum != 0
    expr_R = [1.0 * i / expr_sum for i in expr]
    expr_Log = []
    for i in expr_R:
        if i != 0:
            expr_Log.append(i*ln(i)/ln(2))
        else:
            expr_Log.append(i)
    shannon = -1 * sum(expr_Log)
    return shannon
#-------------------------------

#def sortExprMatrix(exprMatrix):
#    '''
#    exprMatrix = [['id', expr, anno], 
#                  ['id', expr, anno]]
#    '''

#----END sortExprMatrix-------------------------

def output(aDict, matrixD, prefix, annoD):
    '''
    aDict = {'compare': {
                'samp':[]
                'id':[]}   
            }
    matrixD = {'head':[],
             'gene':{'samp':expr, 'samp2':expr2}}
    '''
    for key, valueD in aDict.items():
        output = prefix+key+".results"
        if annoD:
            output_anno = prefix+key+".anno.xls"
            anno_fh = open(output_anno, 'w')
        fh = open(output, 'w')
        sampleL = valueD['samp']
        idL = valueD['id']
        headerL = matrixD['head']
        existL = [samp for samp in headerL if samp in sampleL]
        print >>fh, "%s\t%s\tshannonEntropy" % ('gene', '\t'.join(existL))
        #print >>sys.stderr, matrixD['CG11790']
        if annoD:
            print >>anno_fh, "%s\t%s\tshannonEntropy\t%s" % \
                ('gene', '\t'.join(existL), annoD['head'])
        exprMatrix = []
        for id in idL:
            tmpL = [id]
            exprL = [matrixD[id][samp] for samp in existL]
            tmpL.append('\t'.join(exprL))
            shannon = "%.3f" % computeShannon(exprL)
            tmpL.append(str(shannon))
            if annoD:
                tmpL.append(annoD.get(id, ""))
            #------END one id of each file------------
            exprMatrix.append(tmpL)
        #----------END each item-----------------------
        #exprMatrix = sortExprMatrix(exprMatrix)
        exprMatrix.sort(key=lambda x: x[2])
        for tmpL in exprMatrix:
            print >>fh, '\t'.join(tmpL[:-1])
        #print >>fh, '\n'.join(exprMatrix)
        if annoD:
            print >>anno_fh, \
                '\n'.join(['\t'.join(tmpL) for tmpL in exprMatrix])
        fh.close()
        anno_fh.close()
#---------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    prefix = file.replace("all", "")
    desp = options.filein_desp
    despD = {}
    if desp:
        for line in open(desp):
            value, key = line.split()
            despD[key] = value
    #-----------------------------------------
    expr = options.expr
    samp = options.samp
    sampD = readSamp(samp)
    anno = options.anno
    annoD = readAnno(anno)
    remove = [i.strip() for i in options.remove.split(',')]
    rep  = [i.strip() for i in options.replace.split(',')]
    sep = [i.strip() for i in options.sep.split(',')]
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    for line in fh:
        value, key = line.strip().split('\t')
        oriKey = despD.get(key, key)
        if key not in aDict:
            aDict[key] = {}
            sampL = parseKey(oriKey, sep, remove, rep, sampD)
            aDict[key]['samp'] = sampL
            aDict[key]['id'] = [value]
        else:
            aDict[key]['id'].append(value)
    if debug:
        for key, valueD in aDict.items():
            print >>sys.stderr, "%s\t%s\t%s" %\
                (key, valueD['samp'], valueD['id'][0])
    #-----------------------------------
    matrixD = readMatrix(expr)
    output(aDict, matrixD, prefix, annoD)
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


