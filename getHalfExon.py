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
Functional description

This is used to get junction centered regions. This program get right
half of left exon and left half of right exon.

Input file format:
chr7	52823164	52823749	NR_038166_2.Coding_exon.1   0	-
chr7	52826355	52826562	NR_038166_2.Coding_exon.2   0	-
chr7	52829782	52829892	NR_038166_2.Coding_exon.3   0	-
chr7	52829977	52830147	NR_038166_2.Coding_exon.4   0	-
chr7	52830496	52830546	NR_038166_2.Coding_exon.5   0	-
chr5	31351045	31351129	NM_027855_3.Coding_exon.1   0	+
chr5	31351834	31351953	NM_027855_3.Coding_exon.2   0	+
chr5	31354569	31354641	NM_027855_3.Coding_exon.3   0	+
chr5	31354834	31354906	NM_027855_3.Coding_exon.4   0	+
chr5	31355135	31355257	NM_027855_3.Coding_exon.5   0	+
chr5	31356333	31356431	NM_027855_3.Coding_exon.6   0	+

'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    desc = ""
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="An input file format")
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
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    for line in fh:
        lineL = line.split()
        name  = lineL[3]
        gene, null, num = name.split('.')
        num = int(num)
        if gene not in aDict:
            aDict[gene]= {}
        if num not in aDict[gene]:
            aDict[gene][num] = lineL
        else:
            print >>sys.stderr, "Duplicate num %s %d" % (gene, num)
            sys.exit(1)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    for gene, numD in aDict.items():
        numDKl = numD.keys()
        numDKl.sort()
        lennumDKl = len(numDKl) - 1
        for i in range(lennumDKl):
            num_L = numDKl[i]
            num_R = numDKl[i+1]
            leftL = numD[num_L][:]
            rightL = numD[num_R][:]
            lenLeftExon = (int(leftL[2]) - int(leftL[1])) / 2
            lenRightExon = (int(rightL[2]) - int(rightL[1])) / 2
            strand = leftL[5]
            if strand == '+':
                left_no = '@1'
                right_no = '@2'
            elif strand == '-':
                left_no = '@2'
                right_no = '@1'
            #-------------------------------------------
            leftL[1] = str(int(leftL[2])-lenLeftExon)
            leftL[3] = gene + '.junc.' + str(i+1) + left_no
            print '\t'.join(leftL)
            rightL[2] = str(int(rightL[1])+lenRightExon)
            rightL[3] = gene + '.junc.' + str(i+1) + right_no
            print '\t'.join(rightL)
    #---------------------------------
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



