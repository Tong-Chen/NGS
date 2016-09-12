#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2013, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Functional description:
    Transfer a tow column file to a matrix like displayed below.
Input:
    #Type 1
    pattern name   sequence name
    MA0139.1        A15
    MA0139.1        A25
    MA0139.1        A7
    MA0139.1        A14
    MA0139.1        A13
    MA0139.1        A12
    MA0139.1        A9
    #Type 2
    pattern name   sequence name
    MA0139.1        A15,A1,A3,A4,A7
    MA0139.2        A15,A8,A2,A5,A7

Output
pattern A1  A2  A3  ... A25
MA0139.1    0   0   0   ... 1
UP00021_1   0   0   0   ... 0

'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

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
        metavar="FILEIN", help="A two columns file. \
The first column will be used as the rowname in output file. \
The second column will be used as the colname in output file. \
<-> represents STDIN.")
    parser.add_option("-H", "--header", dest="header",
        default=1, help="Number of header lines to skip. Default 1.")
    parser.add_option("-l", "--absent-present", dest="ap",
        default=0, help="Default only show absent present \
information. If any non-zero number given, will output the \
count of each element in first column grouped by the second \
column.")
    parser.add_option("-s", "--sep-second-col", dest="sep",
        help="Separtor for second column if multiple values contained \
like example input file type 2.")
    parser.add_option("-S", "--skip-blank-item-in-second-col",
        dest="skip_b", default=0, 
        help="Normally blank items in second column after splitting \
represents wrong duplication of separtors. Default <0> meaning \
not skip for backtrack. Normally this parameter should be set \
to <1> to trigger the skipping.")
    parser.add_option("-I", "--full-colname", dest="fullcol",
        help="This file contains the full list of items in \
the second column of file given to -i. This is optional, only \
used when you want to fill 0 to all unmentioned cells.")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    header = int(options.header)
    ap = options.ap
    sep = options.sep
    file2 = options.fullcol
    debug = options.debug
    verbose = options.debug
    skip_b = int(options.skip_b)
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    rownameS = set()
    colnameS = set()

    for line in fh:
        if header:
            header -= 1
            continue
        rowname, colname = line.split()
        if sep:
            colname = colname.split(sep)
            rownameS.add(rowname)
            for each_colname in colname:
                if skip_b:
                    each_colname = each_colname.strip()
                    if not each_colname:
                        continue
                if each_colname not in aDict:
                    aDict[each_colname] = {}
                if rowname not in aDict[each_colname]:
                    aDict[each_colname][rowname] = 0
                aDict[each_colname][rowname] += 1
                colnameS.add(each_colname)
            #----------------------------------------
        else:
            if colname not in aDict:
                aDict[colname] = {}
            if rowname not in aDict[colname]:
                aDict[colname][rowname] = 0
            aDict[colname][rowname] += 1
            #--------------------------------
            rownameS.add(rowname)
            colnameS.add(colname)
            
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    rownameS = list(rownameS)
    rownameS.sort()
    colnameS = list(colnameS)
    colnameS.sort()
    if file2:
        colnameS = [line.strip() for line in open(file2)]
    print "%s\t%s" % ('pattern', '\t'.join(colnameS))
    for rowname in rownameS:
        tmpL = [rowname]
        for colname in colnameS:
            if colname in aDict:
                value = aDict[colname].get(rowname, 0)
                if value != 0 and ap == 0:
                    value = 1
            else:
                value = 0
            #--------------------------
            tmpL.append(str(value))
        #------------------------------------
        print "\t".join(tmpL)
    #-------------------------------------------
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



