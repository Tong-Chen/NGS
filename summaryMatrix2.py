#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import unicode_literals
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
    This is designed to summary matrix statistics.


Normal matrix

ID  sampA   sampB   sampC   ...
id1 1000    1200    1300    ...
id2 1000    1200    1300    ...
id3 1000    1200    1300    ...
id4 1000    1200    1300    ...
id5 1000    1200    1300    ...


Melted format (only two columns allowed)
value   type
1000    sampA
1000    sampA
1000    sampA
1000    sampA
1200    sampB
1200    sampB
1200    sampB
1200    sampB
1200    sampB
1300    sampC
1300    sampC
1300    sampC
1300    sampC
1300    sampC

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import pandas as pd
import re
from tools import targetOlder, targetNewer
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
        metavar="FILEIN", help="A matrix file")
    parser.add_option("-m", "--melted", dest="melted",
        default=False, metavar="MELTED", 
        action="store_true", help="Specify if file in melted format.")
    parser.add_option("-g", "--groupBy-col", dest="groupBy",
        default="type", 
        help="Column name used for grouping (only used for melted format). Default <type>.")
    parser.add_option("-v", "--value-col", dest="value",
        default="value", 
        help="Column name for value column. Default <value>.")
    parser.add_option("-H", "--header", dest="header",
        default=0, type='int', help="An integer to specify the header line. Default <O> represents using first line as header line. Give <-1> to specify no header line." )
    parser.add_option("-r", "--index-col", dest="index_col",
        default='0', help="Column to use as the row labels of the DataFrame. Default <0> meaning the first column as row index. If a sequence like <0,1,2> is given, a MultiIndex will be used. Supply <None> to turn off index_col.")
    parser.add_option("-c", "--usecols", dest="usecols",
        default="None", help="Return a subset of the columns. All elements in this array must either \
be positional (i.e. integer indices into the document columns) or strings \
that correspond to column names provided either by the user in `names` or \
inferred from the document header row(s). For example,  a valid `usecols` \
parameter would be <0,1,2> or <foo,bar,baz>. Using this parameter \
results in much faster parsing time and lower memory usage. \
*******One important note, <index_col> should be also included in <usecols>. \
Default <None> to read in all columns.")
    parser.add_option("-s", "--summarycols", dest="sumcols",
        default="None", help="Specify columns (,separated strings) used for summary analysis. Default all non-index columns.")
    parser.add_option("-x", "--index-name", dest="index_name",
        default="sum", help="Name for the first column in output file. Only works for single column index. Default <sum>. Any string would be OK.")
    parser.add_option("-R", "--remove-all-zero", dest="rm_all0",
        default=False, action="store_true", 
        help="Remove rows containing only zero. Default False.")
    parser.add_option("-o", "--outputfile", dest="output",
        help="Full name for output file. ")
    parser.add_option("-V", "--verbose", dest="verbose",
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
    output = options.output
    if targetNewer(output, file):
        return 0
    melted = options.melted
    groupBy = options.groupBy
    value = options.value
    header = options.header
    if melted:
        assert header == 0, "Header line needed for melted format"
    if header == -1:
        header = 'None'
    index_col = options.index_col.strip()
    if index_col != "None":
        index_col = [int(i) for i in re.split(r'[, ]*', index_col)]
        if len(index_col) == 1:
            index_col = index_col[0]
    else:
        index_col = None
    #-----------------------------------------------------
    usecols = options.usecols
    index_name = options.index_name
    rm_all0 = options.rm_all0
    #print >>sys.stderr, usecols
    if usecols != "None":
        usecols = usecols.split(',')
        usecols_2 = [int(i) for i in usecols if i.isdigit()]
        if len(usecols_2) == len(usecols):
            usecols = usecols_2
    else:
        usecols = None
    #print >>sys.stderr, usecols
    #---------------------------------------------------------
    sumcols = options.sumcols
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    matrix = pd.read_table(file, header=header, index_col=index_col, usecols=usecols)

    #matrix = pd.concat(matrixL, axis=1)
    if debug:
        print >>sys.stderr, matrix.head()
    if isinstance(index_col, int):
        matrix.index.name = index_name
        if debug:
            print >>sys.stderr, matrix.head()
    if sumcols != "None":
        sumcols =  re.split(r'[, ]*', sumcols)
        matrix = matrix.loc[:, sumcols]
    if rm_all0:
        ## Remove all ZEROs
        matrix = matrix.loc[(matrix>0).any(axis=1)]
        if debug:
            print >>sys.stderr, matrix.head()
    if melted:
        matrixGrped = matrix.groupby(groupBy, as_index=True)
        matrixGrped_describ = matrixGrped.describe()
        matrixGrped_describ = matrixGrped_describ.unstack()
        describe = matrixGrped_describ[value].T
    else:
        describe = matrix.describe()
    describe.index.name = index_name
    describe.to_csv(output, sep=b"\t")
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


