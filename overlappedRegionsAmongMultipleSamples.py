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
'''
Functionla description

This is designed to compute the percantage of peaks overlapped among
samples.

It will return a matrix like this.
    a   b   c   d
a   1   0.5 0.7 0.9
b   0.6 1   0.8 0.9
c   0.6 0.4 1   0.5
d   0.7 0.7 0.6 1

Each number in the matrix represents # of regions overlapped between X
and Y relative to total numbe of Y-axis.

This depends on the command <intersectBed> from Bedtools for region
data in bed format and local script <diffComMultiple.py> for one
column name data. 
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#import subprocess

def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    desc = ""
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="Input file must be separated by \
comma only, like 'file1.bed,file2.bed,file3.bed'")
    parser.add_option("-l", "--label", dest="label",
        metavar="FILEIN", help="Unique string to represent \
input files, like 'file1,file2,file3'. Label must be consistent with \
file name. Default using filenmae as labels. If the length of label \
less than the length of filenames, the filenames will be trunctated to \
from tail.")
    parser.add_option("-t", "--type", dest="type",
        metavar="Region/Name", 
        default='Region', help="This indicates the type of data. \
Region represents bed type data, which needs intersectBed to \
compute the numbe of overlapped regions. \
Name represents data with one column or all colums can be \
taken as one thing. This depends on <diffComMultiple.py>. \
Default Region.")
    parser.add_option("-r", "--return-type", dest="returnType",
        metavar="Number/Percentage", 
        default="Percentage", help="The type of data one want to get in \
the output file. Default Percentage.")
    parser.add_option("-n", "--na", dest="NA",
        metavar="1/0", 
        default="1", help="Us NA value for diagonal (left-bottom to \
right-top). When -r is percentage,  default NA value is used. When -r \
is number,  no NA will be used.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    assert options.type in ['Region', 'Name'], \
        "Wrong parameter for -t or --type"
    assert options.returnType in ['Number', 'Percentage'], \
        "Wrong parameter for -r or --return-type"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    fileL = file.split(',')
    lenFileL = len(fileL)
    assert lenFileL > 1, "At least two files are needed"
    label = options.label
    if label == None:
        labeL = fileL
    else:
        labeL = label.split(',')
        lenLabel = len(labeL)
        if lenLabel < lenFileL:
            lenFileL = lenLabel
            print >>sys.stderr, "Truncated file list for less labels \
are given. This may be an error if you do not mean to do this."
            fileL = fileL[:lenFileL]
    labelDict = dict(zip(fileL, labeL))
    verbose = options.verbose
    debug = options.debug
    type = options.type
    returnType = options.returnType
    na = options.NA
    if returnType == 'Percentage' and not na:
        na = 1
    opDict = {'Region':'intersectBed', 'Name':'diffComMultiple.py'}
    aDict = {}
    #-----------------------------------
    for i in fileL:
        for j in fileL:
            if type == 'Region':
                cmd = ["intersectBed -a",i,'-b',j, "-u | wc -l"]
            else:
                cmd = ["diffComMultiple.py",i,j, " | wc -l"]
            #--------------------------------------------------
            aDict[(i,j)] = int(os.popen(' '.join(cmd)).read())
    #-----------------------------------
    print "Name\t%s" % '\t'.join(labeL)
    for i in fileL:
        tmpL = [labelDict[i]]
        total = aDict[(i, i)]
        for j in fileL:
            if returnType == 'Number':
                value = aDict[(i,j)]
            else:
                value = aDict[(i,j)]*1.0/total
                if na and i==j:
                    value = 'NA'
            tmpL.append(str(value))
        #----------------------------------
        print '\t'.join(tmpL)
    #-------------------------------------
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



