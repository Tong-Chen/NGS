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

Input file is 1-based.

Notice:
    bdg: zero-based, half-open.
    wig: 1-based, no strand specific usually
    bed: zero-based, half-open.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
import re

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


def output(chr, adict, col):
    '''
    adict={chr:{pos:value}}
    '''
#-----------------------------

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "\n\t%prog -i file\n\tzcat file.gz | %prog -i -"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="File output by bismark_methylation_extractor_split.py. Gzipped file should be sent in by STDIN.")
    parser.add_option("-o", "--output-prefix", dest="op_prefix",
        help="A string represents output file. The program will add '.wig' suffix and transform to 'bigwig'.")
    parser.add_option("-g", "--genome-prefix", dest="genome_size",
        help="A two-column file with first column as chromosome names and second column as chromosome sizes.")
    parser.add_option("-c", "--coverage-col", dest="coverage_col",
        default=5, type='int', 
        help="Specify the column containing coverage data. Default 5 meaning the fifth column.")
    parser.add_option("-t", "--coverage-threshold", dest="coverage_thres",
        default=10, type='int', 
        help="Specify the coverage threshold for confident sites. Default 10.")
    parser.add_option("-m", "--methylRate-col", dest="methylrate_col",
        default=6, type='int', 
        help="Specify the column containing methylation rate. Default 6 meaning the sixth column.")
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
    verbose = options.verbose
    coverage_thres = options.coverage_thres
    coverage_col   = options.coverage_col - 1
    methylrate_col = options.methylrate_col - 1
    op_prefix      = options.op_prefix
    genome_size    = options.genome_size
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    header = 1
    op_wig = op_prefix + '.wig'
    op_wig_fh = open(op_wig, 'w')
    op_bigwig = op_prefix + '.bw'
    '''
    chrD = {'chr1':1, 'chr2':1}
    '''
    chrD = {}
    ignoreD = {}
    main_chr = re.compile(r'chr[0-9XYMxym]')
    for line in fh:
        if header:
            header -= 1
            continue
        lineL = line.split()
        coverage = int(lineL[coverage_col])
        if coverage < coverage_thres:
            continue

        chr = lineL[1]
        if not main_chr.match(chr):
            if chr not in ignoreD:
                print >>sys.stderr, '{} ignored'.format(chr)
                ignoreD[chr] = 1
            continue
        pos = lineL[2]
        if chr not in chrD:
            chrD[chr] = 1
            print >>op_wig_fh, 'variableStep chrom=%s' % chr

        value = float(lineL[methylrate_col])/100
        #----------------------------------
        print >>op_wig_fh, "{}\t{}".format(pos, value)
    #---------------------------------
    #-------------END reading file----------
    op_wig_fh.close()
    if genome_size:
        cmdL = ['wigToBigWig -clip', op_wig, genome_size, op_bigwig]
        os.system(' '.join(cmdL))
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


