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

This is designed to merge the output of ceas with --dump. 
Files like macs.ceas_dump_TTS.txt, macs.ceas_dump_TSS.txt,
macs.ceas_dump_gene.txt should be given and will be processed.
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
    parser.add_option("-a", "--tss-file", dest="tss",
        metavar="ceas_dump_TSS.txt", help="The output file of ceas \
ended with ceas_dump_TSS.txt.")
    parser.add_option("-b", "--gene-file", dest="gene",
        metavar="ceas_dump_gene.txt", help="The output file of ceas \
ended with ceas_dump_gene.txt.")
    parser.add_option("-c", "--tts-file", dest="tts",
        metavar="ceas_dump_TTS.txt", help="The output file of ceas \
ended with ceas_dump_TTS.txt.")
    parser.add_option("-e", "--extend", dest="extend",
        metavar="rel-sidt for ceas", help="The value given to \
--rel-dist in ceas. In ceas this value defaulted as 3000. \
However, no default value here. Any value given would be transferred \
to integer.")
    parser.add_option("-r", "--pf-res", dest="res",
        metavar="pf-res for ceas", help="The value given to \
--pf-res in ceas. In ceas this value defaulted as 50. \
However, no default value here. Any value given would be transferred \
to integer.")
    parser.add_option("-g", "--num-bins", dest="gene_bin",
        metavar="gene_bin", help="The number of bins \
for each gene. In ceas this value defaulted as 61 (with Tss and Tts \
site involved) when --pf-res is 50. Otherwise you can use \
3000/<--pf-res> +1 to get the value. \
However, no default value here. Any value given would be transferred \
to integer.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    #assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    tss = options.tss
    gene = options.gene
    tts = options.tts
    extend = int(options.extend)
    res = int(options.res)
    num_gene = int(options.gene_bin)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    #----get the head-------------
    headL = []
    num_tss_tts = extend / res #get the number of bins for Tss and Tts
    total_bin = num_tss_tts + num_gene + num_tss_tts
    start = -1 * extend
    for i in range(total_bin):
        headL.append(str(start))
        start += res
    #-----------------------------
    blank_tss_tts = '\t'.join(['0' for i in range(num_tss_tts)])
    blank_gene = '\t'.join(['0' for i in range(num_gene)])
    #-----read tss pat-----------------
    tss_d = {}
    for line in open(tss):
        lineL = line.strip().split("\t")
        name='-'.join(lineL[:4])
        valueL = lineL[5].strip(',').split(',')
        #print valueL
        if len(valueL) > 1:
            tss_d[name] = '\t'.join(valueL[:num_tss_tts])
        else:
            tss_d[name] = blank_tss_tts
    #-------------END reading file----------
    #-----read gene pat-----------------
    gene_d = {}
    for line in open(gene):
        lineL = line.strip().split("\t")
        name='-'.join(lineL[:4])
        valueL = lineL[5].strip(',').split(',')
        if len(valueL) > 1:
            gene_d[name] = '\t'.join(valueL[:num_gene])
        else:
            gene_d[name] = blank_gene
    #-------------END reading file----------
    #-----read tts pat-----------------
    tts_d = {}
    for line in open(tts):
        lineL = line.strip().split("\t")
        name='-'.join(lineL[:4])
        valueL = lineL[5].strip(',').split(',')
        if len(valueL) > 1:
            tts_d[name] = '\t'.join(valueL[num_tss_tts+1:])
        else:
            tts_d[name] = blank_tss_tts

    #-------------END reading file----------
    #----prepare output-----
    len_tss_d = len(tss_d)
    len_gene_d = len(gene_d)
    len_tts_d = len(tts_d)
    assert len_tss_d == len_gene_d and len_tss_d == len_tts_d
    print "Gene\t%s" % '\t'.join(headL)
    geneNl = gene_d.keys()
    geneNl.sort()
    for geneN in geneNl:
        print "%s\t%s\t%s\t%s" % \
            (geneN, tss_d[geneN], gene_d[geneN], tts_d[geneN])
    #--------------------------
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



