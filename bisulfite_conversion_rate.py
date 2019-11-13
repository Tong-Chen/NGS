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
    This is designed to compute conversion rate of bisulfite sequencing.

    When treated with bisulfite, unmodified C (in single strand) 
    will be deaminated and converted into U(T), leaving mC and 5hmC as C.

    Normally ummethylated Lambda DNA will be added to record the conversion rate,  which is the number of unmethylated C divides by all detected C.
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
        metavar="FILEIN", help="Single nucleotide methylation rate file with header line. <-> represents STDIN.")
    parser.add_option("-l", "--sampleLabel", dest="label",
        metavar="sampleLabel", help="Name for the sample.")
    parser.add_option("-m", "--metRate-col", dest="metrate_col",
        type="int", default=8, help="Specify the column containing MetRate to determine the status of cytosines. Default 8 meaning the 7th column.")
    parser.add_option("-t", "--total-cov-col", dest="total_cov_col",
        type="int", default=5, help="Specify the column containing Total coverage to determine the coverage of cytosines. Default 5 meaning the 4th column.")
    parser.add_option("-r", "--met_rate", dest="met_rate",
        default="0", help="Any methylation rate smaller than given float value (default 0, accept values from 0-1) will be treated as unmethylated sites. Deplated.")
    parser.add_option("-c", "--chr-lambda", dest="chr_lambda",
        default="chrLambda", help="Name for lambda DNA, default <chrLambda>. Only lines beginning with this will be counted.")
    parser.add_option("-H", "--header-line", dest="header",
        type="int", default=1, help="Number of heder lines to skip. Default <1>.")
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
    label = options.label
    metrate_col = options.metrate_col-1
    total_cov_col = options.total_cov_col-1

    met_rate = options.met_rate
    met_rate_float = float(met_rate)
    chr_lambda = options.chr_lambda
    header = options.header
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    cntD = {}
    total = 0
    unmethy = 0
    for line in fh:
        if header:
            header -= 1
            continue
        if line.startswith(chr_lambda):
            lineL = line.strip().split('\t')
            met_rate = float(lineL[metrate_col])
            total_cov = float(lineL[total_cov_col])
            total += total_cov
            unmethy += total_cov * (1-met_rate)
            #if met_rate == '0':
            #    if real_met_rate == '0':
            #        unmethy += 1
            #elif float(real_met_rate) < met_rate_float:
            #    unmethy += 1
    print "Conversion_rate\tSample"
    if total:
        con_rate = "{:.2f}".format(unmethy / total * 100)
    else:
        con_rate = "NA"
    print "{}\t{}".format(con_rate, label)
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


