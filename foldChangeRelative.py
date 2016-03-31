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
    This is designed to compute log2 fold change for a data matrix. 

For example we have a data file like:

KO  Tr  Gene    T0  T1  T2  T3
KO1 Tr1 Gene1   0   0   0   1
KO2 Tr2 Gene2   1   1   4   10
KO3 Tr3 Gene3   1   1   4   10
KO4 Tr4 Gene4   1   0   4   10

foldChangeRelative.py -i data -c 4 

foldChangeRelative.py -i data -c 4,5

foldChangeRelative.py -i data -c 4 -l 1

foldChangeRelative.py -i data -c 4 -a 1

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from bs4 import BeautifulSoup
from math import log

#reload(sys)
#sys.setdefaultencoding('utf8')

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
        metavar="FILEIN", help="Data matrix file with head lines.")
    parser.add_option("-c", "--control-column", dest="cc",
        metavar="CONTROL-COLUMN", help="A number like 4 \
representing the forth column as control. A comma separated list \
like 4,5 representing using the forth column and fifth column \
as contaol sample respectively. All columns after will divide \
(when -l is 0, default) or subtract (when -l is 1) the contaol \
column(s).")
    parser.add_option("-l", "--log2-or-not", dest="log2_in",
        default=0, help="Please specify 1 if numbers \
are log2 transformed.")
    parser.add_option("-a", "--added-value", dest="add_value",
        default=1, help="Add specified value to each number to \
avoid dividing by zero.")
    parser.add_option("-L", "--log2-output", dest="log2_out",
        default=1, help="Default fold change will be \
log2 transformed in output to clearly showing the up and down \
times. Specify 0 to get raw fold change.")
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
    control_col = [int(i)-1 for i in options.cc.split(',')]
    log2_in = int(options.log2_in)
    add_value = float(options.add_value)
    log2_out = int(options.log2_out)
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    head = 1
    for line in fh:
        lineL = line.strip().split('\t')
        if head:
            for control_index in control_col:
                control_name = lineL[control_index]
                for other in lineL[control_index+1:]:
                    lineL.append(other+"_vs_"+control_name)
            head -= 1
            print '\t'.join(lineL)
            continue
        #--------------------------------------------------
        for control_index in control_col:
            control_value = float(lineL[control_index])+add_value
            for other_value in lineL[control_index+1:]:
                other_value = float(other_value) + add_value
                if log2_in:
                    diff = other_value-control_value
                else:
                    diff = other_value / control_value
                if log2_out and (diff != 0):
                    diff = log(diff) / log(2)
                lineL.append(str(diff))       
        print '\t'.join(lineL)
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


