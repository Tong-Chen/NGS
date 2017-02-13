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
    This is designed to paste multiple files using pandas.
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
        metavar="FILEIN", help="A list of files with first column as index and first row as header. Files can be separated by <,> or < >.")
    parser.add_option("-l", "--label", dest="label",
        metavar="LABEL", help="A list of names with format and order as FILEIN. If <--rename_cols> is None,  this parameter is not needed.")
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
    parser.add_option("-n", "--rename_cols", dest="rename_cols",
        default="infer", help="Rename columns before merging. If only one column is \
extracted (excluding row index), column name will be renamed with <file_label>; \
If more than one columns are extracted, column names will be renamed with <file_label_original_name>. \
If a list of names are given when extracting multiple columns, \
column names will be renamed with <file_label_given_name>.\
Supply <None> to disable rename.")
    parser.add_option("-m", "--method", dest="method",
        default="outer", help="Join method,  default <outer>, accept <inner>, <left>,  <right>.")
    parser.add_option("-N", "--na", dest="na",
        help="A value given to substitute NA value. Default no substitute.")
    parser.add_option("-R", "--remove-all-zero", dest="rm_all0",
        default=False, action="store_true", 
        help="Remove rows containing only zero. Default False.")
    parser.add_option("-I", "--interger-output", dest="integer",
        default=False, action="store_true", 
        help="When joining matrixes, missed value will be replaced by NaN firtst and other values will be transferred to float. If you still want interger output as the input. Turn on this parameter. Default False.")
    parser.add_option("-o", "--outputfile", dest="output",
        help="Name for output file. If <.gz> suffix is given, compressed file will be output.")
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
    fileL = re.split(r'[, ]*', file.strip())
    if options.rename_cols != "None":
        label = options.label
        labelL = re.split(r'[, ]*', label.strip())
    else:
        labelL = fileL
    header = options.header
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
    na = options.na
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
    rename_cols = options.rename_cols
    if rename_cols not in ["infer", "None"]:
        rename_cols = re.split(r'[, ]*', rename_cols.strip())
    method = options.method
    output = options.output
    integer = options.integer
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    matrixL = []
    #print fileL
    count = 0
    for file, label in zip(fileL, labelL):
        if debug:
            print >>sys.stderr, file
            print >>sys.stderr, header
            print >>sys.stderr, index_col
            print >>sys.stderr, usecols

        matrix = pd.read_table(file, header=header, index_col=index_col, usecols=usecols)
        if isinstance(rename_cols, list):
            column_name = rename_cols
            column_name = [label+'_'+i for i in column_name]
            matrix.columns = column_name
        elif rename_cols != "None":
            column_name = matrix.columns
            if len(column_name) == 1:
                column_name = [label]
            else:
                column_name = matrix.columns
                column_name = [label+'_'+i for i in column_name]
            matrix.columns = column_name
            #print >>sys.stderr, column_name
            if debug:
                print >>sys.stderr, matrix.head()
        matrixL.append(matrix)
        #matrix.to_csv("test."+str(count), sep=b'\t')
        count += 1

    matrix = matrixL[0]
    matrix = matrix.join(matrixL[1:], how=method)
    #matrix = pd.concat(matrixL, axis=1)
    if debug:
        print >>sys.stderr, matrix.head()
    if index_col != None:
        matrix.index.name = "ID"
    if na:
        try:
            na = int(na)
        except ValueError:
            na = float(na)
        except ValueError:
            na = na
        #print >>sys.stderr, na
        matrix = matrix.fillna(na)
        if debug:
            print >>sys.stderr, matrix.head()
    if integer:
        matrix = matrix.astype(int)
    if rm_all0:
        ## Remove all ZEROs
        matrix = matrix.loc[(matrix>0).any(axis=1)]
        if debug:
            print >>sys.stderr, matrix.head()
    describe = matrix.describe()
    if output.endswith('.gz'):
        matrix.to_csv(output, sep=b"\t", compression='gzip')
        describe.to_csv(output.replace('.gz', '.sta'), sep=b"\t")
    else:
        matrix.to_csv(output, sep=b"\t")
        describe.to_csv(output+'.sta', sep=b"\t")
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


