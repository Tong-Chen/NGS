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

This requires all files have same number of rows and same labels in
the first column.

Test passed compared with paste --- 20131111

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
Any number of ',' and <space> is allowed" )
    parser.add_option("-l", "--file-label", dest="file_label",
        metavar="FILE_LABEL", help="The labels should have the same \
format and order as filenames given to <-i>. \
The order matters since a file will \
containling labels by row will be output.")
    parser.add_option("-o", "--output", dest="output",
        metavar="OUTPUT", help="The name of output file containing \
pasted files. The order of rows in output file is the same as first file. \
Also a file named <OUTPUT>.label will be output as described in -l.")
    parser.add_option("-e", "--num-extra-col", dest="num_extra_col",
        default=0, help="The number of extra columns to add. Default 0 means no adding.")
    parser.add_option("-v", "--val-extra-col", dest="val_extra_col",
        default='0', help="The values in extra columns.")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    assert options.file_label != None, "A filelabel needed for -l"
    assert options.output != None, "Output file needed for -o"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    fileL = re.split('[, ]*', options.filein.strip())
    file_labelL = re.split('[, ]*', options.file_label.strip())
    #file_labelL = options.file_label.split(',')
    output = options.output
    num_extra_col = int(options.num_extra_col)
    if num_extra_col:
        val_extra_col = '\t'.join([options.val_extra_col] *
            num_extra_col)
    #--------------------------------------------

    fh = open(output+".label", 'w')
    print >>fh, '\n'.join(file_labelL)
    fh.close()
    debug = options.debug
    #-----------------------------------
    aDict = {}
    labelL = []
    i = 0
    lenFileL = len(fileL)
    for file in fileL:
        if debug:
            print >>sys.stderr, file
        aDict[file] = {}
        i += 1
        for line in open(file):
            lineL = line.strip().split("\t", 1)
            key = lineL[0]
            if i == 1:
                labelL.append(key)
            if debug:
                print >>sys.stderr, lineL
            value = lineL[1]
            if i < lenFileL and num_extra_col:
                value = '\t'.join([value, val_extra_col])
            if key not in aDict[file]:
                aDict[file][key] = value
            else:
                print "Duplicate keys %s in file %s" % (key, file)
                sys.exit(1)
        #---------------------------------------------------
    #-------------END reading file----------
    fh = open(output, 'w')
    for key in labelL:
        print >>fh, "%s\t%s" \
            % (key, '\t'.join([aDict[file][key] for file in fileL]))
    fh.close()
    #-----------end close fh-----------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()



