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
    This is designed to normalize a matrix.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import pandas as pd
import re
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
        metavar="FILEIN", help="A numerical matrix file. With first line as header line and first row as index row.")
    parser.add_option("-T", "--transpose-data", dest="transpose",
        default=False, action="store_true", 
        help="Transpose data before normalize.")
    parser.add_option("-r", "--index-col", dest="index_col",
        default='0', help="Column to use as the row labels of the DataFrame (Default 0 represents the first column). If a sequence like <0,1,2> is given, a MultiIndex will be used. Supply <None> to turn off index_col.")
    parser.add_option("-c", "--usecols", dest="usecols",
        default="None", help="Return a subset of the columns. All elements in this array must either \
be positional (i.e. integer indices into the document columns) or strings \
that correspond to column names provided either by the user in `names` or \
inferred from the document header row(s). For example,  a valid `usecols` \
parameter would be <0,1,2> or <foo,bar,baz>. Using this parameter \
results in much faster parsing time and lower memory usage. \
*******One important note, <index_col> should be also included in <usecols>. \
Default <None> to read in all columns.")
    parser.add_option("-n", "--norm_factor", dest="norm_factor",
        help="A file containing normalize factor for each column in data matrix. \
The row name of norm_factor should be the same as column name of data matrix (order does not matter). Specially a string <sum> will use the sum value of each column as \
normalize factor. A string <StandardScaler> or <MinMaxScaler> will do column-wise normalization as indicated.")
    parser.add_option("-o", "--outputfile", dest="output",
        help="[Lowercase o] Output file prefix. Default input file with suffix <.xls, .tsv, .csv, .txt removed if exists>. A suffix <xls.gz> will be added.")
    parser.add_option("-u", "--uncompressed", dest="uncompressed",
        default=False, action="store_true", help="Default output gzipped file. Specify to output uncompressed result.")
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
    transpose = options.transpose
    header = 0
    index_col = options.index_col.strip()
    if index_col != "None":
        index_col = [int(i) for i in re.split(r'[, ]*', index_col)]
    else:
        index_col = None
    #-----------------------------------------------------
    norm_factor = options.norm_factor
    usecols = options.usecols
    if usecols != "None":
        usecols = usecols.split(',')
        usecols_2 = [int(i) for i in usecols if i.isdigit()]
        if len(usecols_2) == len(usecols):
            usecols = usecols_2
    else:
        usecols = None
    #---------------------------------------------------------
    uncompressed = options.uncompressed
    if not options.output:
        suffix = ['xls', 'tsv', 'csv', 'txt']
        if file.find('.') != -1:
            file_name, suf = file.rsplit('.', 1)
            if suf not in suffix:
                file_name = file
        else:
            file_name = file
        output = file_name
    else:
        output = options.output
    assert output, "-o/--outputfile should be supplied."
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    matrix = pd.read_table(file, header=header, index_col=index_col, usecols=usecols)
    if transpose:
        matrix = matrix.T
    matrix.index.name = "ID"
    if norm_factor in ['StandardScaler', 'MinMaxScaler']:
        from sklearn import preprocessing
        matrix = matrix.astype('float')
        x = matrix.values
        if norm_factor == 'StandardScaler':
            scaler = preprocessing.StandardScaler()
        elif norm_factor == 'MinMaxScaler':
            scaler = preprocessing.MinMaxScaler()
        x_scaled = scaler.fit_transform(x)
        matrix2 = pd.DataFrame(x_scaled)
        matrix2.columns = matrix.columns
        matrix2.index   = matrix.index
        matrix = matrix2
    else:
        if norm_factor == 'sum':
            norm_variable = matrix.sum()
        else:
            normD = pd.read_table(norm_factor, header=None, index_col=0)
            norm_variable = normD[1]
        matrix = matrix / norm_variable
    #--------------------------------------------------------
    if uncompressed:
        full = output + '.xls'
        matrix.to_csv(full, sep=b"\t")
    else:
        full = output + '.xls.gz'
        matrix.to_csv(full, sep=b"\t", compression='gzip')
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


