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

This is designed to paste multiple files and remove labels of
non-first file and confirm all lines have same label as listed in the
first column.

Test passed compared with paste --- 20131111
'''

import sys
import os
import re
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
        metavar="FILEIN", help="Multiple files can be separated with \
<,> or space(< >) or both, like <'file1,file2,file3'>, \
<'file1 file2 file3'> or <'file1' 'file2', 'file3'>." )
    parser.add_option("-l", "--file-label", dest="file_label",
        metavar="FILE_LABEL", help="The labels should have the same \
format and order as filenames given to <-i>. \
Nomally the label will be added as another columns of original file \
to indicate the source of corrsonding rows. \
However this parameter is not necessary.")
    parser.add_option("-c", "--col", dest="column",
        default=0, help="Select the columns you want to output. \
Normally the first column is label column, and will only be output \
once. Default 0 means output all columns of all files. \
Accept '2,3,4' or '2-5' or '2,3-6,7' to indicate the columns wanted. \
'2,3-6,7' means column 2-7 will all be output.")
    parser.add_option("-r", "--row-list", dest="row_list",
        help="Blank separated lists (normally a subset of the first \
column in each file, other items not in given files are also \
accepted. However their value wolud be 0.) \
represent the items you want to output and their order.")
    parser.add_option("-o", "--output", dest="output",
        metavar="OUTPUT", help="The name of output file containing \
pasted files. The line order of output file is the same as first file.")
    parser.add_option("-H", "--HEADER", dest="header",
        default=0, help="When <1> is given, skip header lines for all files \
except the first one.")
    parser.add_option("-n", "--NEW_HEADER", dest="new_header",
        help="Set header line explicitly by giving ',' or ' ' \
separated strings like 'variable,value,Sample'. Default using header lines of \
the first file.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    fileL = re.split('[, ]*', options.filein.strip())
    #fileL = options.filein.split(',')
    col = options.column
    if col == 0:
        colL = []
    else:
        colL = []
        for i in col.split(','):
            if i.find('-') == -1:
                colL.append(int(i)-1)
            else:
                range1 = i.split('-')
                for j in range(int(range1[0]), int(range1[1])+1):
                    colL.append(j-1)
            #------------------------------------
        #--------------------------------------
    #---------------------------------------------
    output = options.output
    header = int(options.header)
    if options.file_label:
        #file_labelL = options.file_label.split(',')
        file_labelL = re.split('[, ]*', options.file_label.strip())
    #--------------------------------
    verbose = options.verbose
    debug = options.debug
    new_header = options.new_header
    if new_header:
        new_header = '\t'.join(re.split('[, ]*', new_header.strip()))
    rowlist = options.row_list
    if rowlist:
        rowlist = re.split('[, ]*', rowlist.strip())
    #-----------------------------------
    aDict = {}
    labelL = []
    i = 0
    headerL =[]
    for file in fileL:
        aDict[file] = {}
        i += 1
        tmpSkip = header
        for line in open(file):
            if tmpSkip:
                tmpSkip -= 1
                if i == 1:
                    headerL.append(line.strip())
                continue
            #---------------------------
            specifyCol = len(colL)
            if specifyCol == 0:
                lineL = line.strip().split("\t", 1)
                key = lineL[0]
                if i == 1:
                    labelL.append(key)
                value = lineL[1]
            else:
                lineL = line.strip().split("\t")
                key = lineL[0]
                if i == 1:
                    labelL.append(key)
                    #---------------------------------------
                valueL = [lineL[j] for j in colL]
                value = '\t'.join(valueL)
            #--------------------------------------
            if key not in aDict[file]:
                aDict[file][key] = value
            else:
                print "Duplicate keys %s in file %s" % (key, file)
                sys.exit(1)
        #---------------------------------------------------
    #-------------END reading file----------
    #print aDict
    #--------------------------
    fh = open(output, 'w')
    if new_header:
        print >>fh, "%s" % new_header
    else:
        if header and headerL:
            if options.file_label:
                print >>fh, "%s\tSample" % '\t'.join(headerL)
            else:
                print >>fh, "%s" % '\t'.join(headerL)
        #-----------------------------
    #----------------------------------
    if rowlist:
        labelL = rowlist
    
    i = -1
    for file in fileL:
        i += 1
        for label in labelL:
            if options.file_label:
                print >>fh, "%s\t%s\t%s" % (label,
                        aDict[file].get(label, '0'), file_labelL[i])
            else:
                print >>fh, "%s\t%s" % (label, aDict[file].get(label, '0'))
    #--------------------------------------
    #for key in labelL:
    #    print >>fh, "%s\t%s" \
    #        % (key, '\t'.join([aDict[file].get(key, '0') for file in fileL]))
    fh.close()
    #-----------end close fh-----------
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



