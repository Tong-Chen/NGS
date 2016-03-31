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

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

debug = 0

def fprint(content):
    print json_dumps(content,indent=1)

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i cuffcompare.gtf -f cuff -g gname -t trname >output"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="GTF file")
    parser.add_option("-f", "--gtf-type", dest="gtf_type",
        metavar="FILEIN", help="Specific the type of gtf file, \
currently accept, <cuff> representing GTFs generated using cuff* \
tools; <dexseq> meansing GTFs generated for DEXSeq.")
    parser.add_option("-g", "--gname", dest="gname",
        help="Two columns file with the first column as \
<gene_id> in GTF file and second column as <substitute_gene_id>. \
Normally this can be generated using <parseCuffmerge.merged_gtf.py>.")
    parser.add_option("-t", "--trname", dest="trname",
        help="Two columns file with the first column as \
<transcript_id> in GTF file and second column as \
<substitute_transcript_id>. Normally this can be generated using \
<parseCuffmerge.merged_gtf.py>.")
    parser.add_option("-s", "--strand", dest="strand",
        default='no', help="Default <no> indicating no \
considering of strand information. If <yes> is given, \
misambiguous transcripts will be removed.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
#def readNameTransfer(file):
#    aDict = {}

#------------------------------------------
def parseGTF(gtf, gtf_type, gnameD, trnameD, strand):
    if gtf_type == "cuff":
        for line in open(gtf):
            lineL = line.split('\t')
            if strand == 'yes' and lineL[6] == ".":
                continue
            nameAttrL = lineL[-1].strip(';\n').split(';')
            nameD = dict([item.strip('" ').split(' "') for item in nameAttrL])
            gene_id = nameD['gene_id']
            newGene_id = gnameD.get(gene_id, gene_id)
            tr_id = nameD['transcript_id']
            newTr_id = trnameD.get(tr_id, tr_id)
            line = line.replace(gene_id, newGene_id)
            line = line.replace(tr_id, newTr_id)
            print line,
    elif gtf_type == 'dexseq':
        for line in open(gtf):
            lineL = line.split('\t')
            if strand == 'yes' and lineL[6] == ".":
                continue
            nameAttrL = lineL[-1].strip(';\n').split(';')
            nameD = dict([item.strip('" ').split(' "') for item in nameAttrL])
            gene_id = nameD['gene_id']
            gene_idL = gene_id.split('+')
            newGene_id = '+'.join([gnameD.get(i,i) for i in gene_idL])
            line = line.replace(gene_id, newGene_id)
            tr_id = nameD.get('transcripts', ' ')
            if tr_id != ' ':
                newTr_id = '+'.join([trnameD.get(i,i) for i in tr_id.split('+')])
                line = line.replace(tr_id, newTr_id)
            print line,
    else:
        print >>sys.stderr, "Unsupported GTF type"

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    gtf = options.filein
    gtf_type = options.gtf_type
    gname = options.gname
    gnameD = dict([line.strip().split('\t') for line in open(gname)])
    trname = options.trname
    trnameD = dict([line.strip().split('\t') for line in open(trname)])
    verbose = options.verbose
    global debug
    debug = options.debug
    strand = options.strand
    #-----------------------------------
    parseGTF(gtf, gtf_type, gnameD, trnameD, strand)
    #-----------end close fh-----------
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
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
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


