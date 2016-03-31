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
    This is designed to get the expression of genes.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
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
        metavar="FILEIN", help="Input matrix")
    parser.add_option("-f", "--function-file", dest="func_file",
        metavar="FUNC_FILE", help="One column file containing \
specific functional terms one want to get.")
    parser.add_option("-l", "--label-file", dest="label_file",
        metavar="LABEL_FILE", help="One column file containing \
specific IDs one want to label.")
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
    func_file = options.func_file
    if func_file:
        funcL = [line.strip().upper() for line in open(func_file)]
        funcNameL = [line.strip().replace(' ', '_') for line in open(func_file)]
    else:
        funcL = []
    lenfuncL = len(funcL)
    label_file = options.label_file
    if label_file:
        labelL = [line.strip().upper() for line in open(label_file)]
    else:
        labelL = []
    lenlabelL = len(labelL)
    #-------------------------------------------------
    verbose = options.verbose
    global dbug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    header = 1
    for line in fh:
        line = line.strip()
        if header:
            if funcL:
                line += '\tFunc'
            if labelL:
                line += '\tLabel'
            print line
            header -= 1
            continue
        #--------------------------
        lineU = line.upper()
        output = 0
        if funcL:
            findFirst = 0
            for i in range(lenfuncL):
                func = funcL[i]
                if lineU.find(func) != -1:
                    if not findFirst:
                        line += '\t'+funcNameL[i]
                    else:
                        line += ';'+funcNameL[i]
                    findFirst = 1
                    output = 2
                    #break
        #--------------------------------------------
        if labelL:
            find_label = 0
            for label in labelL:
                if lineU.find(label) != -1:
                    line += '\t*'
                    find_label = 1
                    break
            if not find_label:
                line += '\t-'
            if funcL and output == 2:
                output = 1
            if not funcL:
                output = 1
        #------------------------------------
        if output:
            print line
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


