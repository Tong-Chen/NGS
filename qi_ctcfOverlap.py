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
        metavar="FILEIN", help="")
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
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    for line in fh:
        lineL = line.strip().split('\t')
        insulator = lineL[3]
        label = lineL[4]
        if insulator not in aDict:
            aDict[insulator] = {}
        if label not in aDict[insulator]:
            aDict[insulator][label] = [lineL]
        else:
            aDict[insulator][label].append(lineL)
        #---------------------------------------------------
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    for insulator, innerD in aDict.items():
        #print insulator, innerD
        if "CTCF" in innerD:
            ctcfL = innerD["CTCF"]
        else:
            ctcfL = [['chr_test_ct', '-10', '-1', 'CT_CTCF', 'CTCF', 'MA0139.1']]
        for key, valueL in innerD.items():
            if key != 'CTCF':
                for element in valueL:
                    output = 1
                    elestart = int(element[1])
                    eleend   = int(element[2])
                    eleSet = set()
                    for i in range(elestart, eleend):
                        eleSet.add(i)
                    lenEleSet = len(eleSet)
                    for ctcf in ctcfL:
                        ctcfstart = int(ctcf[1])
                        ctcfend   = int(ctcf[2])
                        ctcfSet = set()
                        for i in range(ctcfstart, ctcfend):
                            ctcfSet.add(i)
                        #---------------------------------------
                        lenCtcfSet = len(ctcfSet)
                        overlap = len(ctcfSet.intersection(eleSet))
                        if overlap * 1.0 / lenEleSet > 0.5:
                            output = 0
                            break
                    #---------------------------------------
                    if output:
                        print '\t'.join(element)
                    else:
                        if verbose:
                            print ctcfL
                            print element
            #-------------------------------
        #-------------------------------------------

    #---------------------------------------
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



