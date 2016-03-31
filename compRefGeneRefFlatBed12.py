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

    This script is used to compare the refGene.gtf and refFlat.gtf to
    generate standard GTF file.

    refGene.gtf:
    chr1    hg19_refGene    start_codon     67000042        67000044    0.000000        +       .       gene_id " NM_032291" ;  transcript_id " NM_032291" ; 
    chr1    hg19_refGene    CDS     67000042        67000051    0.000000        +       0       gene_id " NM_032291" ;  transcript_id " NM_032291" ; 
    chr1    hg19_refGene    exon    66999825        67000051    0.000000        +       .       gene_id " NM_032291" ;  transcript_id " NM_032291" ;
    refFlat.gtf:
    chr1    hg19_refGene    start_codon     67000042        67000044    0.000000        +       .       gene_id " SGIP1" ;  transcript_id " SGIP1" ; 
    chr1    hg19_refGene    CDS     67000042        67000051    0.000000        +       0       gene_id " SGIP1" ;  transcript_id " SGIP1" ; 
    chr1    hg19_refGene    exon    66999825        67000051    0.000000        +       .       gene_id " SGIP1" ;  transcript_id " SGIP1" ;

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
    parser.add_option("-i", "--refGene", dest="filein",
        metavar="refGene.bed12", help="")
    parser.add_option("-f", "--refFlat", dest="refFlat",
        metavar="refFlat.bed12", help="")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readFile(fh):
    aDict = {}
    for line in fh:
        lineL = line.strip().split('\t')
        value = lineL[3]
        lineL.remove(value)
        keyL = '\t'.join(lineL)
        aDict[keyL] = value
    #-----------------------------------------------
    return aDict
#---------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    refFlat = options.refFlat
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    refGeneD = readFile(fh)
    refFlatD = readFile(open(refFlat))
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    #for key1, valueL1 in refGeneD.items():
    #    for key2, valueL2 in refFlatD.items():
    #        if valueL1 == valueL2:
    #            print "%s\t%s" % (key1, key2)
    #            refFlatD.pop(key2)
    #            break
    #----------------------------------------------
    for key1, valueL1 in refGeneD.items():
        valueL2 = refFlatD[key1]
        print "%s\t%s" % (valueL1, valueL2)


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



