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
    This is designed to transform a multiple column to a matrix with
    one column as row names, one column as column names and one column
    as values.
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
        metavar="FILEIN", help="A file contains at least three same \
length columns separated with tab.")
    parser.add_option("-r", "--row-name-column", dest="row",
        default=1, help="Specify which column used as row names in \
outputted matrix. Default 1 means the first column.")
    parser.add_option("-c", "--column-name-column", dest="col",
        default=2, help="Specify which column used as column names in \
outputted matrix. Default 2 means the second column.")
    parser.add_option("-V", "--value-column", dest="value",
        default=3, help="Specify which column used as values in \
outputted matrix. Default 3 means the third column.")
    parser.add_option("-H", "--header-line", dest="header",
        default=1, help="The number of header lines one want to skip. \
Default 1.")
    parser.add_option("-l", "--full-list-row-name", dest="full_row",
        help="A file contains full list of row names.")
    parser.add_option("-m", "--missing-value", dest="miss_value",
        default="NA", help="Specifying missing values. Default <NA>.")
    parser.add_option("-D", "--duplicate_names", dest="dpn",
        default=1, help="A number to represent the number of duplicate \
names in name column. Default 1 means no duplicate names allowed. \
A specified integer indicating the number of duplicate names.")
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
    row = int(options.row) - 1
    col = int(options.col) - 1
    value = int(options.value) - 1
    miss_value = options.miss_value
    header = int(options.header)
    full_row = options.full_row
    dpn = int(options.dpn)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    sampleSet = set()
    aDict = {}
    for line in fh:
        if header:
            header -= 1
            continue
        #-------------------------------
        lineL  = line.strip().split("\t")
        gene   = lineL[row]
        sample = lineL[col]
        expr   = lineL[value]
        sampleSet.add(sample)
        if gene not in aDict:
            aDict[gene] = {}
        if dpn > 1:
            if sample not in aDict[gene]:
                aDict[gene][sample] = [expr]
            else:
                aDict[gene][sample].append(expr)
        elif dpn == 1:
            if sample not in aDict[gene]:
                aDict[gene][sample] = expr
            else:
                print >>sys.stderr, "Duplicate sample for gene", gene, sample
    #-------------END reading file----------
    #-------------Output--------------------
    if dpn > 1:
        sampleSet = list(sampleSet) 
        sampleSet.sort()
        sampleSet2 = sampleSet * dpn
        sampleSet2.sort()
        print "Name\t%s" % '\t'.join(sampleSet2)
        for gene, sampleD in aDict.items():
            tmpL = ['\t'.join(sampleD.get(sample, [miss_value]*dpn)) for sample in sampleSet]
            print "%s\t%s" % (gene, '\t'.join(tmpL))
    #---------------------------------------------------------------------
    elif dpn == 1:
        sampleSet = list(sampleSet)
        sampleSet.sort()
        print "Name\t%s" % '\t'.join(sampleSet)
        if full_row:
            geneL = [line.strip() for line in open(full_row)]
        else:
            geneL = aDict.keys()
        for gene in geneL:
            sampleD = aDict.get(gene, {})
        #for gene, sampleD in aDict.items():
            tmpL = [sampleD.get(sample, miss_value) for sample in sampleSet]
            print "%s\t%s" % (gene, '\t'.join(tmpL))
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



