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
    This is designed to transfer pandas summary file for sp_boxplot_precomputed_values.sh to get boxplot.

Input format

	Au_113_F	Au_113_M	Au_113_C
count	927076.0	927076.0	927076.0
mean	0.776247698062	0.779987818688	0.77837770253
std	0.272107045425	0.27213324813	0.258101695504
min	0.0	0.0	0.0
25%	0.571428571429	0.571428571429	0.55737704918
50%	0.909090909091	0.923076923077	0.932432432432
75%	1.0	1.0	1.0
max	1.0	1.0	1.0

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import pandas as pd
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
        metavar="FILEIN", help="Format as specified above")
    parser.add_option("-p", "--pheno-file", dest="pheno",
        help="Sample information. Optional")
    parser.add_option("-o", "--output-file", dest="output",
        help="Output file name")
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
    pheno = options.pheno
    verbose = options.verbose
    output = options.output
    global debug
    debug = options.debug
    #-----------------------------------
    data = pd.read_table(file, header=0, index_col=0)
    data.index.name = "ID"
    data = data.T
    renameD = {"min":"minimum", "max":"maximum", "25%":"lower_quantile", 
            "50%":"median", "75%":"upper_quantile"}
    data.rename(columns=renameD,  inplace=True)
    if pheno:
        phenoData = pd.read_table(pheno, header=0, index_col=0)
        dataM = data.join(phenoData, how="left")
    else:
        dataM = data
    dataM.index.name = "ID"
    dataM.to_csv(output, sep=b"\t")
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


