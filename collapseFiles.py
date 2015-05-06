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
Functionla description

This is designed to paste multiple files (no matter the number of
columns) and remove labels of
non-first file and confirm all lines have same label as listed in the
first column of the first file.

There is no need for all files containing same number or order of
labels. The program will colllect all labels to generate the final
label list.

Optionally, you can specify one file containing all labels as the row_label
file to get wanted row order.

'''

import sys
import os
from time import localtime, strftime 
import re
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
        metavar="FILEIN", help="Multiple files can be separated with \
<,> or <space> like <'file1,file2,file3'> or <'file1 file2 , file3'>. \
Any number of ',' and <space> is allowed. The program does not \
assume all files have same number of rows in each file or these rows \
have same order. It will automatically generate a full list of row \
labels (based on the first column of each file) in alphabetical order. \
Or you can supply a full list of row labels to <-m> to specify the \
order of output files.")
    parser.add_option("-l", "--file-label", dest="file_label",
        metavar="FILE_LABEL", help="The labels should have the same \
format and order as filenames given to <-i>. They will appear as \
header line to indicate the source of each column. This is not necessary.")
    parser.add_option("-L", "--row-label", dest="row_label",
        help="A file containing the full list of row names. If this is \
given, the order of rows in output file will be the same as in this \
<row_label> file.")
    parser.add_option("-m", "--missing-value", dest="miss",
        default='0', help="A value used to pad the matrix when not all files \
have same content in the first column. Default 0.")
    parser.add_option("-o", "--output", dest="output",
        metavar="OUTPUT", help="The name of output file containing \
pasted files. Default STDOUT. The order of rows in output file is alphabetically \
sorted or the same as the row_label file when given. ")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A list of filenames needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    fileL = re.split('[, ]*', options.filein.strip())
    if options.file_label:
        file_labelL = re.split('[, ]*', options.file_label.strip())
    else:
        file_labelL = ''
    #----------------------------------------------
    if options.output:
        output = options.output
        out_fh = open(output, 'w')
    else:
        out_fh = sys.stdout
    #----------------------------------------------
    if options.row_label:
        labelL = [i.strip() for i in open(options.row_label)]
        generate_label = 0
    else:
        labelL = set()
        generate_label = 1
    #-----------------------------------
    debug = options.debug
    missing = options.miss
    #-----------------------------------
    aDict = {}
    lenFileL = len(fileL)
    
    i = -1
    for file in fileL:
        if debug:
            print >>sys.stderr, file
        aDict[file] = {}
        i += 1
        header = 1
        for line in open(file):
            lineL = line.strip().split("\t", 1)
            key = lineL[0]
            if generate_label:
                labelL.add(key)
            if debug:
                print >>sys.stderr, lineL
            value = lineL[1]
            #----------Deal with missing value and header lines-----
            if header:
                value_len = value.count('\t') + 1
                aDict[file]['miss'] = '\t'.join([missing] * value_len)
                if file_labelL:
                    aDict[file]['header'] = \
                        '\t'.join([file_labelL[i]] * value_len)
                header -= 1
            #----------Deal with missing value and header lines-----
            if key not in aDict[file]:
                aDict[file][key] = value
            else:
                print "Duplicate keys %s in file %s" % (key, file)
                sys.exit(1)
        #---------------------------------------------------
    #-------------END reading file----------
    if generate_label:
        labelL = list(labelL)
        labelL.sort()
    if file_labelL:
        print >>out_fh, "%s\t%s" \
            % ('Name', '\t'.join([aDict[file].get('header', 'miss_header') for file in fileL]))
    #----------------------------------------------------------------------
    for key in labelL:
        print >>out_fh, "%s\t%s" \
            % (key, '\t'.join([aDict[file].get(key,
                aDict[file]['miss']) for file in fileL]))
    if options.output:
        out_fh.close()
    #-----------end close fh-----------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()



