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
    This is designed to transfer a matrix to three columns file.

Input:
    ID  A   B   C   D
    A   1   2   3   5
    B   2   3   8   9
    C   2   0   0   2

Output:
    A   A   1
    A   B   2
    A   C   3
    A   D   5
    B   A   2
    B   B   3
    B   C   8
    .
    .
    .
    C   D   2

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

def fprint(content):
    """ 
    This is a Google style docs.

    Args:
        param1(str): this is the first param
        param2(int, optional): this is a second param
            
    Returns:
        bool: This is a description of what is returned
            
    Raises:
        KeyError: raises an exception))
    """
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
        metavar="FILEIN", help="A matrix file as described above.")
    parser.add_option("-s", "--symmeric", dest="symmeric",
        default=False, action="store_true", help="If the matrix is symmeric, \
only half information will be output. Default False.")
    parser.add_option("-t", "--threshold", dest="threshold",
        help="Limit output for pairs which value meets given threshold limitation. \
Accepted format like a single numeric value or ',' spearated two values like '-0.8,0.8'.")
    parser.add_option("-c", "--compute-way", dest="comp_way",
        help="For single neumric value, one can give <less> (output pairs which value smaller than threshold), <lessequal> (output pairs which value smaller than or equal threshold), <larger>, <largerequal>. For paired values, one can give <inner> (output pairs which value between given two values), <outer> (output pairs which value at outer sides of given two values), <innerequal>, <outerequal>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def determineThreshold(value, threshold, len_threshold, comp_way):
    if threshold:
        value = float(value)
        if len_threshold == 1:
            threshold = threshold[0]
            assert comp_way in ['less', 'lessequal', 'larger', 'largerequal']
            if comp_way == 'less' and value<threshold:
                return 1
            elif comp_way == 'lessequal' and value <=threshold:
                return 1
            elif comp_way == 'larger' and value > threshold:
                return 1
            elif comp_way == 'largerequal' and value >= threshold:
                return 1
            else:
                return 0
        elif len_threshold == 2:
            threshold.sort()
            assert comp_way in ['inner', 'innerequal', 'outer', 'outerequal']
            if comp_way == 'inner' and threshold[0]<value<threshold[1]:
                return 1
            elif comp_way == 'innerequal' and threshold[0]<=value<=threshold[1]:
                return 1
            elif comp_way == 'outer' and (threshold[0]>value or threshold[1]<value):
                return 1
            elif comp_way == 'outerequal' and (threshold[0]>=value or threshold[1]<=value):
                return 1
            else:
                return 0
        #-----------------------------------------            
        else:
            print >>sys.stderr, "Wrong threshold format"
            sys.exit(1)
    #-----------------------------------------            
    return True
#--------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    symmeric = options.symmeric
    threshold = options.threshold
    len_threshold = 0
    if threshold:
        threshold = [float(i) for i in threshold.split(',')]
        len_threshold = len(threshold)
    comp_way = options.comp_way

    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    header = 1
    saveD = {}
    for line in fh:
        lineL = line.split()
        if header:
            headerL = lineL[1:]
            header -= 1
            continue
        key = lineL[0]
        saveD[key] = 1
        for seckey, value in zip(headerL, lineL[1:]):
            if seckey == key:
                continue
            if not symmeric or seckey not in saveD:
                if debug:
                    print >>sys.stderr, value, threshold, comp_way
                output = determineThreshold(value, threshold, len_threshold, comp_way)
                if debug:
                    print >>sys.stderr, output
                if output:
                    print "\t".join([key, seckey, value])
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


