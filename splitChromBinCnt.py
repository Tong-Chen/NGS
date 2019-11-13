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
    This is designed to split bin size cnt file by chromosomes. Chromosome name is determined by the first column. Those string before the last underline will be treated as chromosome names.

Input file:
    chr1_7959       8
    chr1_7960       6
    chr1_7961       2
    chr1_7962       1
    chr1_7964       2
    chr2_7965       2
    chr2_7966       2
    chr2_7967       2
    chr2_7968       2
    chrUn_GL000218v1_3207   1
    chrUn_GL000218v1_3215   1
    chrUn_GL000218v1_3216   1
    chrUn_GL000218v1_3217   1
    chrUn_GL000218v1_3218   3
    chrUn_GL000218v1_3219   4

'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
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
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="Input file with format specified above")
    parser.add_option("-c", "--key-column", dest="key_col", 
            type='int', default=1, 
            help="The column contaning keys to be split on. Default 1 meaning the first column.")
    parser.add_option("-a", "--all-chrom", dest="all_chrom", 
            help="A file with first column containing all chromosome names. If specified, the chromome not in <filein> will be also generated as an empty file. Optional.")
    parser.add_option("-o", "--output-prefix", dest="op_prefix",
        help="A string to specify prefix of output files. The output file will be named in format like <prefix>.<chr_name>.xls")
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
    key_col = options.key_col-1
    all_chrom = options.all_chrom
    if all_chrom:
        chromL = [i.strip().split('\t')[0] for i in open(all_chrom)]
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
    chromD = {}
    for line in fh:
        lineL = line.split()
        key = lineL[key_col]
        chrom = key.rsplit('_', 1)[0]
        if chrom not in chromD:
            chromD[chrom] = open(op_prefix+'.'+chrom+'.xls', 'w')
        output_fh = chromD[chrom]
        print >>output_fh, line, 
    #-------------END reading file----------
    for chrom, output_fh in chromD.items():
        if all_chrom and chrom in chromL:
            chromL.remove(chrom)
        output_fh.close()
    #-------------------------------
    #print >>sys.stderr, all_chrom
    #print >>sys.stderr, chromL
    if all_chrom and chromL:
        for chrom in chromL:
            os.system("touch {}.{}.xls".format(op_prefix, chrom))
            #print >>sys.stderr, "touch {}.{}.xls".format(op_prefix, chrom)
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


