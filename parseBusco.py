#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
from __future__ import print_function
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to parse BUSCO short_summary result.
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print(desc, file=sys.stderr)
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="Short summary")
    parser.add_option("-o", "--output-file", dest="op_prefix",
        metavar="FILEIN", help="Output file")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    op_prefix = options.op_prefix
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {'C': 'Complete BUSCOs', 'S':'Complete and single-copy BUSCOs', 'D':'Complete and duplicated BUSCOs', 
            'F': 'Fragmented BUSCOs', 'M': 'Missing BUSCOs'}
    numD = {}
    for line in fh:
        if line[0] == '#':
            continue
        if line.find('%') != -1:
            line = line.replace('[', ',')
            line = line.replace(']', '')
            percentD = dict([i.split(':') for i in line.strip().split(',')])
            continue           
        if line[0] == "\t":
            num, str1 = line.strip().split('\t')
            par = str1.find('(')
            if par != -1:
                label = str1[par+1:par+2]
                if label != 'C':
                    numD[label] = num
    #-------------END reading file----------
    if file != '-':
        fh.close()
    with open(op_prefix, 'w') as op_fh:
        print("x\tvariable\tvalue\tset", file=op_fh)
        for key, value in percentD.items():
            value = value.replace('%', '')
            if key not in ['C', 'n']:
                tmpL = '\t'.join(['x', aDict[key], value, 'Percent'])
                print(tmpL, file=op_fh)
        
        for key, value in numD.items():
            if key != 'C':
                tmpL = '\t'.join(['x', aDict[key], value, 'Raw number'])
                print(tmpL, file=op_fh)
    #----close file handle for files-----
    cmd = ['sp_barPlot.sh -f', op_prefix, 
            "-m TRUE -a x -P right -B set -O 1 -X FALSE -k 'free_y' -N TRUE"]
    os.system(' '.join(cmd))
    #-----------end close fh-----------
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------

