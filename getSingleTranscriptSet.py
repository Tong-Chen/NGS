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

Get the single transcript set for each gene. For genes with multiple
transcripts, the longest one will be picked out.

Input file format(only Coding_exon and UTRs are used):

chr7    52823749    52829782    NR_038165_1.Intron.1    0   -
chr7    52829892    52829977    NR_038165_1.Intron.2    0   -
chr7    52830147    52830496    NR_038165_1.Intron.3    0   -
chr7    52823164    52823749    NR_038165_1.Coding_exon.1   0   -
chr7    52829782    52829892    NR_038165_1.Coding_exon.2   0   -
chr7    52829977    52830147    NR_038165_1.Coding_exon.3   0   -
chr7    52830496    52830546    NR_038165_1.Coding_exon.4   0   -
chr7    52823164    52830546    NR_038165_1 0610005C13Rik   -
chr5    31351045        31351129        NM_027855_3.Coding_exon.1 0       +
chr5    31351834        31351953        NM_027855_3.Coding_exon.2 0       +
chr5    31354569        31354641        NM_027855_3.Coding_exon.3 0       +
chr5    31354834        31354906        NM_027855_3.Coding_exon.4 0       +
chr5    31355135        31355257        NM_027855_3.Coding_exon.5 0       +
chr5    31356333        31356431        NM_027855_3.Coding_exon.6 0       +
chr5    31356635        31356740        NM_027855_3.Coding_exon.7 0       +
chr5    31351012        31356996        NM_027855_3     0610007C21Rik +
chr5    31351012        31351045        NM_027855_3.UTR5        0 +
chr5    31356740        31356996        NM_027855_3.UTR3        0 +

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
        metavar="FILEIN", help="Usually \
mm9.gene.refgene.gtf.normal.bed. Standard input is unsuitable. \
The second column in output file contains R,S,M. \
R: representative transcripts. M: multiple transcripts, S: single \
transcript.")
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
    fh = open(file)
    #--------------------------------
    geneD = {}
    trD   = {}
    trToGene = {}
    '''
    geneD = {geneName: {trname:length}, }
    trD   = {trname:length}
    trToGene = {tr:gene}
    '''
    for line in fh:
        lineL = line.split()
        gene  = lineL[4]
        tr    = lineL[3]
        if gene != '0':
            gene = gene + '@' + lineL[0]
            trToGene[tr] = gene
            if gene not in geneD:
                geneD[gene] = {}
            geneD[gene][tr] = 0
            #else:
            #    print >>sys.stderr, \
            #        "Duplicate names %s" % gene
        #---------------------------------------
    fh.close()
    #---------------------------------
    for line in open(file):
        lineL = line.split()
        gene  = lineL[4]
        keyL  = lineL[3].split('.')
        if gene == '0' and (keyL[1].find('UTR') != -1 or \
            keyL[1].find('Coding_exon') != -1):
            tr = keyL[0]
            gene = trToGene[tr]
            geneD[gene][tr] += int(lineL[2]) - int(lineL[1])
        #------------------------------------------
    #-------------END reading file----------
    #----close file handle for files-----
    #----output--------------------------
    for gene, trD in geneD.items():
        if len(trD) > 1:
            type = 'M'
            tmp_trL = trD.keys()
            tmp_trL.sort()
            tr = tmp_trL[0]
            lenv = trD[tr]
            print "%s\t%s\t%s\t%d" % (tr, 'M', gene, lenv)
            for tmp_tr in tmp_trL[1:]:
                tmp_lenv = trD[tmp_tr]
                print "%s\t%s\t%s\t%d" % (tmp_tr, 'M', gene, tmp_lenv)
                if tmp_lenv > lenv:
                    lenv = tmp_lenv
                    tr   = tmp_tr
            #---------------------------------
            print "%s\t%s\t%s\t%d" % (tr, 'R', gene, lenv)
        else:
            type = 'S'
            print "%s\t%s\t%s\t%d" % (trD.keys()[0], type, gene, \
                trD.values()[0])
        #------------------------
    #-----------end -----------
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



