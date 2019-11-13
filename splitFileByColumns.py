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
    Split files by columns.

For example, if we have a file like below
expr    Control Gene1   Gene2   Gene3
a   1   1   1   1
b   1   1   1   1
c   1   1   1   1
d   1   1   1   1

The program will generate a file for each Gene (for example, Gene1, Gene2)

expr    Control Gene1
a   1   1
b   1   1
c   1   1
d   1   1


expr    Control Gene2
a   1   1
b   1   1
c   1   1
d   1   1

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import numpy as np
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
        metavar="FILEIN", help="A matrix file with the first line as header line")
    parser.add_option("-k", "--kept-cols", dest="kept_cols",
        default=1, help="Columns kept for all subfiles. Default 1 representing the first column. Given <1,2> will keep first two columns for all subfiles.")
    parser.add_option("-o", "--output-prefix", dest="output_prefix",
        help="Prefix for output subfiles.")
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
    kept_cols = [int(i)-1 for i in options.kept_cols.split(',')]
    output_prefix = options.output_prefix
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
    all_indexL = []
    aDict = {}
    for line in fh:
        lineL = line.strip().split('\t')
        if header:
            header -= 1
            len_line = len(lineL)
            all_indexL = range(len_line)
            for index in all_indexL:
                aDict[index] = [lineL[index]]
            continue
        #--------------------------------------
        for index in all_indexL:
            aDict[index].append(lineL[index])
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    kept_values = []
    for index in all_indexL:
        if index in kept_cols:
            kept_values.append(aDict[index])
    for index in all_indexL:
        if index in kept_cols:
            continue
        tmp_keptL = kept_values[:]
        valueL = aDict[index]
        tmp_keptL.append(valueL)
        tmp_keptL = np.array(tmp_keptL).T
        output = output_prefix + '.'+valueL[0]+'.xls'
        tmp_keptL.tofile(output, sep="\t")
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


