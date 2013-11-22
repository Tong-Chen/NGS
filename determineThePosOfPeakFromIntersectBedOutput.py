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

This is designed to determine the position of a peak along one
transcript.

Usually the input file is the output from intersectBed with called
peaks assigned to -a and annotation assigned to -b with -wao as
additional parameters.

Sample lines (tab separated,  exon number in the 12th column is not
required. Currently the transcript name info is needed. Only lines
with matched transcript name in the forth and twelveth column are
used.):
chr14 66772168 66772343 NM_001162366_22522__6__U1 1.73269149085 - UTR3 175 chr14 66772093 66772890 NM_001162365_22521.UTR3 0 - 175
chr14 66772168 66772343 NM_001162366_22522__6__U1 1.73269149085 - UTR3 175 chr14 66772093 66772890 NM_001162366_22522.UTR3 0 - 175
chr14 66772168 66772343 NM_001162366_22522__6__U1 1.73269149085 - UTR3 175 chr14 66772093 66772890 NM_172498_22523.UTR3 0 - 175
chr2 18596576 18596636 NM_147778_6766__1__S__1 1.36528630201 + Coding_exon.7-Coding_exon.8-UTR3 146 chr2 18596576 18596636 NM_147778_6766.Coding_exon.7 0 + 60

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
        metavar="FILEIN", help="The output file of intersectBed. \
Detailed infromation would be showed with -h parameter. STDIN as - \
is accepted.")
    parser.add_option("-p", "--peak_name_col", dest="peak_name_col",
        default=4, help="The column number to represent the column \
harboring peak names. Default 4 (1-based).")
    parser.add_option("-s", "--peak_name_sep", dest="peak_name_sep",
        default="__", help="The separtor used to get the real name \
of peak. Only used when one peak has been splitted into two lines. \
Default '__'. 'FALSE' means unsed.")
    parser.add_option("-n", "--peak_subname_col", dest="peak_subname_col",
        default=3, help="Only used when -s is not 'FALSE'. \
The number of elements extracted to represent the \
peak name. Default 3 means the first three elements.")
    parser.add_option("-m", "--match_type_col", dest="match_type_col",
        default=12, help="The column number to represent the column \
harboring matched types. Default 12 (1-based).")
    parser.add_option("-o", "--output_col", dest="output_col",
        default=8, help="The last column you want to output. \
Every column from 1 to this given number will be outputted. \
Default 8 (1-based) means output 1-8 column.")
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
    file = options.filein
    peak_name_col = int(options.peak_name_col)-1
    peak_name_sep = options.peak_name_sep
    peak_subname_col = int(options.peak_subname_col)
    match_type_col = int(options.match_type_col)-1
    output_col = int(options.output_col)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    lineDict = {}
    typeDict = {}
    for line in fh:
        lineL = line.strip().split('\t')
        full_name = lineL[peak_name_col]
        #if peak_same_sep != 'FALSE':
        peak_nameL = full_name.split("__")
        peak_tr = peak_nameL[0]
        peak_name = "__".join(peak_nameL[:3])
        matched_name = lineL[match_type_col].split(".") 
        matched_tr = matched_name[0]
        matched_type = matched_name[1] 
        if peak_tr == matched_tr:
            if peak_name not in lineDict:
                lineDict[full_name] = lineL[:output_col]
            else:
                assert lineDict[full_name] == lineL[:output_col]
            if peak_name not in typeDict:
                typeDict[peak_name] = {}
            if matched_type not in typeDict[peak_name]:
                typeDict[peak_name][matched_type] = int(lineL[-1])
            else:
                typeDict[peak_name][matched_type] += int(lineL[-1])
        #--------excluding peaks mapped to other transcripts-----
    #-------------END reading file----------
    typeCol = 6
    for full_name, lineL in lineDict.items():
        peak_nameL = full_name.split('__')
        peak_name  = '__'.join(peak_nameL[:3])
        matched_typeL = typeDict[peak_name].keys()
        matched_typeL.sort()
        #if len(matched_typeL)==1:
        #    lineL[typeCol] = matched_typeL[0]
        #    print '\t'.join(lineL)
        #else:
        type = matched_typeL[0]
        type_len = typeDict[peak_name][type]
        for tmpType in matched_typeL[1:]:
            tmpType_len = typeDict[peak_name][tmpType]
            if tmpType_len > type_len:
                type_len = tmpType_len
                type = tmpType
        #--------------------------------------------
        all_type = '.'.join(matched_typeL)
        if type != all_type:
            type +=  '@' + all_type
        lineL[typeCol] = type
        print '\t'.join(lineL)
    #----close file handle for files-----
    if file != '-':
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



