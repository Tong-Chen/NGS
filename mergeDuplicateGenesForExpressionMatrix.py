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
    This is designed to merge genes with multiple probes together.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

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
        metavar="FILEIN", help="Gene expression file")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def average(valueL):
    '''
    valueL = [[1, 2, 3], [4, 5, 6], [], ...]
    '''
    len1 = len(valueL[0])
    len2 = len(valueL)
    tmpL = [0 for i in range(len1)]
    for i in range(len1):
        for value in valueL:
            tmpL[i] += value[i]
        tmpL[i] = str(tmpL[i] / float(len2))
    return tmpL
#---------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    header = 1
    aDict  = {}
    for line in fh:
        if header:
            print line,
            header -= 1
            continue
        #------------------------------
        lineL = line.split()
        gene = lineL[0]
        if gene in aDict:
            if len(aDict[gene]) == 1:
                tmpL = aDict[gene][0]
                aDict[gene] = [[float(i) for i in tmpL]]
            #----------------------------------------------
            aDict[gene].append([float(i) for i in lineL[1:]])
        else:
            aDict[gene] = [lineL[1:]]
        #----------------------------------------------------
    #-------------END reading file----------
    for gene, value in aDict.items():
        if len(value) > 1:
            print '%s\t%s' % (gene, '\t'.join(average(value)))
        else:
            print '%s\t%s' % (gene, '\t'.join(value[0]))
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


