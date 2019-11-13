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
    This is designed to compute fold-change ratio for time-series data.

1. Input matrix file

   Normal expression matrix with rows as genes and columns as samples.

   ID  WT_0h_1 WT_0h_2 WT_1h_1 WT_1h_2 WT_2h_1 WT_2h_2 MT_0h MT_1h
   A    5   5   10  10  20  20  6   9   
   B    15   15   10  10  20  20  15   9   
   C    15   15   10  10  20  20  15   9   
   D    15   15   40  60  10  30  0   9   

2. compare_table 
    (1. every sample in first column will be used as control sample to compute fold changes)
    (2. samples in second column will be used for output order for columns of input matrix.)

   WT_0h    WT_1h
   WT_0h    WT_2h
   MT_0h    MT_1h

3. sampleFile (specify group information for each sample)

   sample    Group
   WT_0h_1  WT_0h
   WT_0h_2  WT_0h
   WT_1h_1  WT_1h
   WT_1h_2  WT_1h
   WT_2h_1  WT_2h
   WT_2h_2  WT_2h
   MT_0h    MT_0h
   MT_1h    MT_1h

'''

def readsampleFile(sampleFile):
    '''
    sampleFileD = {'WT_0h_1': 'WT_0h', 'WT_0h_2', 'WT_0h'}
    '''
    header =1
    sampleFileD = {}
    for line in open(sampleFile):
        if header:
            header -= 1
            continue
        #------------------------------
        sample, grp = line.split()
        sampleFileD[sample] = grp
    return sampleFileD
#-------------------------------------
def readCompare_table(compare_table):
    '''
    compareL = ['WT_1h', 'WT_2h', 'MT_1h']
    compareD = {"WT_1h": "WT_0h", "WT_2h":"WT_0h", "MT_1h":"MT_0h"}
    '''
    compareL = []
    compareD = {}
    for line in open(compare_table):
        control, treat = line.split()
        compareL.append(treat)
        assert treat not in compareD, "Duplicate %s in %s" % (treat, compare_table)
        compareD[treat] = control
    return compareL, compareD
#--------------------------------

import sys
import os
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from math import log

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
        metavar="FILEIN", help="")
    parser.add_option("-s", "--sampleFile", dest="sampleFile",
        help="Format as specified above.")
    parser.add_option("-c", "--compare_table", 
        dest="compare_table",
        help="Format as specified above")
    parser.add_option("-a", "--add-pseudo-count", 
        dest="pseudo_count",
        help="Add a pseudo count value if 0 contains.")
    parser.add_option("-l", "--log2", 
        dest="log2", default=False, action="store_true", 
        help="Specify this parameter if data have already been log2 transformed.")
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
    sampleFile = options.sampleFile
    compare_table = options.compare_table
    pseudo_count = float(options.pseudo_count)
    log2 = options.log2
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    sampleFileD = readsampleFile(sampleFile)
    compareL, compareD = readCompare_table(compare_table)
    header = 1
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    pos2sampD = {}
    for line in fh:
        lineL = line.strip().split()
        if header:
            header -= 1
            len_lineL = len(lineL)
            print "{}\t{}".format(lineL[0], '\t'.join(compareL))
            for i in range(2, len_lineL):
                pos2sampD[i] = lineL[i]
            continue
        #-----save value-------------------------
        valueD = {}
        for i in range(2, len_lineL):
            sample = pos2sampD[i]
            group  = sampleFileD[sample]
            if group not in valueD:
                valueD[group] = []
            valueD[group].append(float(lineL[i]))
        #------get average-------------------------
        for group, valueL in valueD.items():
            valueD[group] = sum(valueL)/len(valueL)
        #------Get ratio--------------------------
        ratioL = []
        for grp in compareL:
            control = compareD[grp]
            control_value = valueD[control]
            grp_value     = valueD[grp]
            if log2:
                ratio = grp_value - control_value
            else:
                if control_value == 0:
                    grp_value += pseudo_count
                    control_value += pseudo_count
                ratio = grp_value / control_value
                ratio = log(ratio)/log(2)
            ratioL.append("%.3f" % ratio)
        print "{}\t{}".format(lineL[0], '\t'.join(ratioL))

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


