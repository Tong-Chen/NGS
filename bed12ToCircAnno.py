#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to transfer bed 12 file for < Gene Predictions and RefSeq Genes with Gene Names>
    format for CIRCexplorer2 usages.

    Output file formamt:

    Field   Description
    geneName    Name of gene
    isoformName Name of isoform
    chrom   Reference sequence
    strand  + or - for strand
    txStart Transcription start position
    txEnd   Transcription end position
    cdsStart    Coding region start
    cdsEnd  Coding region end
    exonCount   Number of exons
    exonStarts  Exon start positions
    exonEnds    Exon end positions

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

def fprint(content):
    """ 
    This is a Google style docs.

    Args:
        param1(str): this is the first param
        param2(int, optional): this is a second param
            
    Returns:
        bool: This is a description of what is returned
            
    Raises:
        KeyError: raises an exception))
    """
    print json_dumps(content,indent=1)

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--gtf-file", dest="filein",
        metavar="FILEIN", help="GTF file")
    parser.add_option("-b", "--bed12", dest="bed12",
        help="GTF transfferred BED12 file")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def geneId_trID_pair(gtf):
    count = 0
    geneId_trID_pairD = {}
    for line in open(gtf):
        lineL = line.split("\t")
        if lineL[2] == 'transcript':
            count += 1
            annoCol = lineL[8]
            annoColL = line.strip(';').split(';')
            for item in annoColL:
                if item.find("gene_id") != -1:
                    gene = item.split('"')[1]
                elif item.find("transcript_id") != -1:
                    transcript = item.split('"')[1]
                #try:
                #    key_t, value_t = item.strip('"').split('"')
                #except ValueError:
                #    print >>sys.stderr, lineL
                #    print >>sys.stderr, item
                #    sys.exit(1)
                #key_t = key_t.strip()
                #value_t = value_t.strip()
                #if key_t == "gene_id":
                #    gene = value_t
                #elif key_t == "transcript_id":
                #    transcript = value_t
            #-------------------------
            if transcript:
                geneId_trID_pairD[transcript] = gene
            else:
                print >>sys.stderr, "Please check yout GTF file. \
                    No 'transcript_id' in the third column found."
                sys.exit(1)
            transcript = gene = ""
    #--------------------------------
    if not count:
        print >>sys.stderr, "Please check yout GTF file. No 'transcript' in the third column found."
        print >>sys.stderr, "We will try to use the slow version"
        geneId_trID_pairD = geneId_trID_pair_slow(gtf)
        #sys.exit(1)
    return geneId_trID_pairD
#-------------------------------------------------

def geneId_trID_pair_slow(gtf):
    count = 0
    geneId_trID_pairD = {}
    for line in open(gtf):
        lineL = line.split("\t")
        count += 1
        annoCol = lineL[8]
        annoColL = line.strip(';').split(';')
        for item in annoColL:
            if item.find("gene_id") != -1:
                gene = item.split('"')[1]
            elif item.find("transcript_id") != -1:
                transcript = item.split('"')[1]
            #try:
            #    key_t, value_t = item.strip('"').split('"')
            #except ValueError:
            #    print >>sys.stderr, lineL
            #    print >>sys.stderr, item
            #    sys.exit(1)
            #key_t = key_t.strip()
            #value_t = value_t.strip()
            #if key_t == "gene_id":
            #    gene = value_t
            #elif key_t == "transcript_id":
            #    transcript = value_t
        #-------------------------
        if transcript:
            geneId_trID_pairD[transcript] = gene
        else:
            print >>sys.stderr, "Please check yout GTF file. \
                No 'transcript_id' in the third column found."
            sys.exit(1)
        transcript = gene = ""
    #--------------------------------
    if not count:
        print >>sys.stderr, "Please check yout GTF file. No 'transcript' in the third column found."
    return geneId_trID_pairD
#-------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    gtf = options.filein
    bed12 = options.bed12
    if not bed12:
        bed12 = gtf + '.bed12'
        if not os.path.exists(bed12):
            if os.system("gtf2bed12.sh -f "+gtf):
                print >>sys, stderr, "Wrong to get bed12"
                sys.exit(1)
    #------------------------------------
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    geneId_trID_pairD = geneId_trID_pair(gtf)
    
    for line in open(bed12):
        lineL = line.split()
        tr = lineL[3]
        start = int(lineL[1])
        exonCnt = int(lineL[-3])
        try:
            exonSizel = [int(i) for i in lineL[-2].strip(',').split(',')]
        except ValueError:
            print >>sys.stderr, lineL[-2]
            sys.exit(1)
        exonPosl  = [int(i) for i in lineL[-1].strip(',').split(',')]
        exonStartl = [start+i for i in exonPosl]
        exonEndl   = [exonStartl[i]+exonSizel[i] for i in range(exonCnt)]
        
        exonStart = ','.join([str(i) for i in exonStartl])+','
        exonEnd = ','.join([str(i) for i in exonEndl])+','

        gene = geneId_trID_pairD.get(tr, tr)
        chr_name = lineL[0]
        newtmpL = [gene, tr, chr_name, lineL[5], lineL[1], lineL[2], lineL[6], 
                lineL[7], lineL[-3], exonStart, exonEnd]
        print '\t'.join(newtmpL)
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


