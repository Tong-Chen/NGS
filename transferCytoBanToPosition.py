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
    This is designed to transfer cytoband postion to chromosome
    postion.


'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import re
#from multiprocessing.dummy import Pool as ThreadPool

def fprint(content):
    print >>sys.stderr, json_dumps(content,indent=1)

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
        metavar="FILEIN", help="The file contains cytoband info \
to transfer.")
    parser.add_option("-c", "--column", dest="cyto_col",
        help="The column containing cytoband information. \
1-based number.")
    parser.add_option("-C", "--cytoband", dest="cytoband",
        help="The cytoband file downloaded from UCSC. \
The file must be sorted first by chromosome and then \
by the position from small to large.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readCytoband(file):
    cytoDict = {}
    #cytoDict = {chr1:{'p36.33': [0, 2300000]}}
    mDict = {}
    #mDict = {chr1: [[min, p], [max, q]]}
    for line in open(file):
        lineL = line.split('\t')
        chr = lineL[0]
        start = int(lineL[1])
        end = int(lineL[2])
        cyto = lineL[3]
        if chr not in mDict:
            mDict[chr] = [[start, cyto[0]], [end, cyto[0]]]
        else:
            if start < mDict[chr][0][0]:
                mDict[chr][0][1] = cyto[0]
            if end > mDict[chr][1][0]:
                mDict[chr][1][1] = cyto[0]
        #---------------------------------------------
        if chr not in cytoDict:
            cytoDict[chr] = {}
        assert cyto not in cytoDict[chr], \
            ' '.join(["Suplicate", cyto, "for", chr, 'in', cytoDict])
        cytoDict[chr][cyto] = [start, end]
        cyto = cyto[:-1].rstrip('.')
        while len(cyto) > 0:
            if cyto in cytoDict[chr]:
                if start < cytoDict[chr][cyto][0]:
                    cytoDict[chr][cyto][0] = start
                if end > cytoDict[chr][cyto][1]:
                    cytoDict[chr][cyto][1] = end
            else:
                cytoDict[chr][cyto] = [start, end]
            cyto = cyto[:-1].rstrip('.')
    #---------END reading-----------------
    for chr, valueLL in mDict.items():
        for valueL in valueLL:
            cyto = valueL[1]+'ter'
            assert cyto not in cytoDict[chr], \
                ' '.join(["Suplicate", cyto, "for", chr, 'in', cytoDict])
            cytoDict[chr][cyto] = [valueL[0]]
    #-------pter and qter------------------
    #fprint(cytoDict)
    return cytoDict
#-----------END readCytoband---------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    cyto_col = int(options.cyto_col) - 1
    cytoband_file = options.cytoband
    verbose = options.verbose
    debug = options.debug

    reg = re.compile('([0-9XYxy]*)([pq][0-9\.ter]*)([pq][0-9ter\.]*)*')

    #-----------------------------------
    cytoDict = readCytoband(cytoband_file)
    #--------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    for line in fh:
        lineL = line.strip().split('\t')
        cyto = lineL[cyto_col].replace('-', '')
        cytoM = reg.match(cyto)
        if cytoM:
            cytoM = cytoM.groups()
            #print cytoM
            unfound = 0
            chr = 'chr'+cytoM[0].upper()
            position = []
            for sub_cyto in cytoM[1:]:
                if not sub_cyto:
                    continue
                if sub_cyto not in cytoDict[chr]:
                    unfound = 1
                    break
                position.extend(cytoDict[chr][sub_cyto])
            if unfound:
                print "%s\t-\t-\t%s" % (chr, line.strip())
                continue
            position.sort()
            print "%s\t%d\t%d\t%s" % (chr, position[0], position[-1],
                    line.strip())
        else:
            print >>sys.stderr, line,
            print >>sys.stderr, "Unracognized cyoband %s" % cyto
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


