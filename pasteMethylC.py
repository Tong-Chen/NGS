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
    This is designed to paste multiple files using pandas specifally for methylation profile.

Input file format

#Chr    Pos     Ref     Chain   Total   Met     UnMet   MetRate Ref_context     Type
chr1    10497   C       +       2       2       0       1       CGG     CpG
chr1    10502   C       +       2       0       2       0       CTG     CHG
chr1    10506   C       +       2       0       2       0       CCT     CHH
chr1    10507   C       +       2       0       2       0       CTG     CHG
chr1    10517   C       +       2       0       2       0       CTG     CHG
chr1    10522   C       +       2       0       2       0       CTC     CHH

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
        metavar="LABEL", help="A list of names with format and order as FILEIN.")
    parser.add_option("-H", "--header", dest="header",
        default=0, type='int', help="An integer to specify the header line. Default <O> represents using first line as header line. Give <-1> to specify no header line." )
    parser.add_option("-r", "--index-col", dest="index_col",
        default='0,1,2,3', help="Column to use as the row labels of the DataFrame. If a sequence like <0,1,2> is given, a MultiIndex will be used. Supply <None> to turn off index_col. Default <0,1,2,3>.")
    parser.add_option("-c", "--usecols", dest="usecols",
        default="0,1,2,3,7,9", help="Return a subset of the columns. All elements in this array must either \
be positional (i.e. integer indices into the document columns) or strings \
that correspond to column names provided either by the user in `names` or \
inferred from the document header row(s). For example,  a valid `usecols` \
parameter would be <0,1,2> or <foo,bar,baz>. Using this parameter \
results in much faster parsing time and lower memory usage. \
*******One important note, <index_col> should be also included in <usecols>******.\
default=<0,1,2,3,7,9>. Supply <None> to readin all columns.")
    parser.add_option("-n", "--rename_cols", dest="rename_cols",
        default="infer", help="Rename columns before merging. If only one column is \
extracted (excluding row index), column name will be renamed with <file_label>; \
If more than one columns are extracted, column names will be renamed with <file_label_original_name>. \
If a list of names are given when extracting multiple columns, \
column names will be renamed with <file_label_given_name>. \
Supply <None> to disable rename.")
    parser.add_option("-m", "--method", dest="method",
        default="outer", help="Join method,  default <outer>, accept <inner>, <left>,  <right>.")
    parser.add_option("-o", "--outputprefix", dest="output",
        help="Name for output prefix.")
    parser.add_option("-C", "--type-col", dest="type_col",
        default="Type", help="Specify the name of column containing methylated C context information like <CpG, CHH,  CHG>. Default <Type> column.")
    parser.add_option("-t", "--type", dest="type",
        default="CpG", help="Default only output <CpG> sites. Accept <all> to output all types os sites.")
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
    label = options.label
    labelL = re.split(r'[, ]*', label.strip())
    header = options.header
    if header == -1:
        header = 'None'
    index_col = options.index_col.strip()
    if index_col != "None":
        index_col = [int(i) for i in re.split(r'[, ]*', index_col)]
    else:
        index_col = None
    #-----------------------------------------------------
    usecols = options.usecols
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
    #-------------------------------------------------------
    type_col = options.type_col
    type = options.type
    if type == 'all':
        typeL = ["CpG", "CHH", "CHG"]
    else:
        typeL = [type]
    method = options.method
    prefix = options.output
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    matrixD = {}
    for type in typeL:
        matrixD[type] = []
    #print fileL
    count = 0
    for file, label in zip(fileL, labelL):
        print >>sys.stderr, file 
        #print >>sys.stderr, header
        #print >>sys.stderr, index_col
        #print >>sys.stderr, usecols

        matrix = pd.read_table(file, header=header, index_col=index_col, usecols=usecols)
        for type in typeL:
            #print >>sys.stderr, type_col
            #print >>sys.stderr, list(matrix.columns)
            keep_col = list(matrix.columns)
            keep_col.remove(type_col)
            #print >>sys.stderr, keep_col
            matrix_tmp = matrix.ix[matrix[type_col]==type, keep_col]

            if isinstance(rename_cols, str) and rename_cols != "None":
                column_names = matrix_tmp.columns
                if len(column_names) > 1:
                    column_name = [label+'_'+i for i in column_names]
                else:
                    column_name = [label]
                matrix_tmp.columns = column_name
            elif isinstance(rename_cols, list):
                column_name = [i+'_'+label for i in rename_cols]
                matrix_tmp.columns = column_name

            matrixD[type].append(matrix_tmp) 
        #matrix.to_csv("test."+str(count), sep=b'\t')
        count += 1

    #matrix = matrixL[0]
    #matrix = matrix.join(matrixL[1:], how=method)
    for type in typeL:
        output = prefix+'.'+type+'.xls.gz'
        matrix_tmpD = matrixD[type]
        matrix_tmp = matrix_tmpD[0].join(matrix_tmpD[1:], how=method)
        matrix_tmp.to_csv(output, sep=b"\t", compression='gzip')
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


