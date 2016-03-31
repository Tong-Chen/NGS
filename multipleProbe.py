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
    This is used to deal with microarray data when one gene
    corresponds to multiple probes.
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
        metavar="FILEIN", help="Gene expression matrix with first \
column representing gene names and first row as header lines.")
    parser.add_option("-o", "--operation", dest="op",
        default='average', help="The way to deal with \
duplicate expressions, default average meaning get the average \
expression of all probes. Accept <max> to keep the probe with the \
biggest difference, and <min> to keep the probe with the \
smalllest difference.")
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
    op   = options.op
    if op != "average":
        print >>sys.stderr, "Currently unsuppored"
        sys.exit(1)
    verbose = options.verbose
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
            print line, 
            header -= 1
            continue
        #--------------------
        lineL = line.split()
        key = lineL[0]
        if key not in aDict:
            aDict[key] = []
        aDict[key].append(lineL[1:])
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    for key, valueLL in aDict.items():
        if len(valueLL) < 2:
            print "%s\t%s" % (key, '\t'.join(valueLL[0]))
        else:
            if op == "average":
                len1 = len(valueLL[0])
                num = len(valueLL)
                newValueL = []
                for i in range(len1):
                    sum = 0.0
                    for valueL in valueLL:
                        sum += float(valueL[i])
                    average = sum / num
                    newValueL.append(str(average))
                #---------------------------------
                print "%s\t%s" % (key, '\t'.join(newValueL))
            #---------END average-----------------------------
            elif op == "max":
                pass
            elif op == "min":
                pass

        #-----------------------------------------------

    #----------------------------------
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


