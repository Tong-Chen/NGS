#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
import sys
import os
from time import localtime, strftime
from re import compile #parseFastqc
#-------global variable------------------------
timeformat          = "%Y-%m-%d %H:%M:%S"
totalReads          = 0
mappedReads         = 0
unmappedReads       = 0
uniqueMappedReads   = 0
duplicateReads      = 0
seed                = 0
qualityEncode       = ''
miniSeqForSeed      = 35 #the minimum reads length of having a seed
seqQuality          = 20 #the minimum standard for a quantified reads 
#-------global variable------------------------
def parseFastqc(fastqc_data):
    '''
    This function parse the output of fastqc, it will return the
    postulated seed, and will promote a line to ask if this seed need
    to be resetted.
    '''
    print "[parseFastqc] Begins %s" % strftime(timeformat, localtime())
    global totalReads
    seedC = 0 #the flag of setting seed 
    for line in open(fastqc_data):
        if line.startswith('Encoding'):
            qualityEncode = line.srip().split('\t')[1]
            if qualityEncode.find('Sanger') != -1:
                qualityEncode = 'sanger'
            else:
                qualityEncode = 'illumina'
            print "[quality encode] detected by fastqc is %s." %\
                qualityEncode
        elif line.startswith('Total Sequences'):
            innerTotalReads =  int(line.split()[1])
            if totalReads:
                try:
                    assert(totalReads == innerTotalReads)
                except AssertionError:
                    print >>sys.stderr, "Unequal total Reads"
            else:
                totalReads = innerTotalReads
            #------------------------------------------
            print "[totalReads] by fastqc is %d." % innerTotalReads
        elif line.startswith('Sequence length'):
            seqlen = int(line.split()[1])
            print "[reads length] by fastqc is %d. \
                miniSeqForSeed is %d." % (seqlen, miniSeqForSeed)
            if seqlen > miniSeqForSeed:
                seedC = 1
                numStart = compile('r[0-9]{2,}')
        elif seedC and numStart.match(line):
            lineL = line.split('\t',3)
            mean = int(lineL[1])
            median = int(lineL[2])
            if mean < seqQuality and \
                median < seqQuality:
                seed = int(lineL[0])
                print "[seed] The chosed seed is %d, with median\
                    %d, mean %d, stardard quality %d " % \
                    (seed, median, mean, seqQuality)
                input = raw_input("Are you satisfied by seed %d?\n \
                If not, please input one integer. \
                If yes, just enter! \n>>>" % seed)
                if input:
                    seed = int(input)
                    assert seed > miniSeqForSeed
                    print "[seed] User set seed is %d" % seed
                #------------------------------------------
                break
        elif line.startswith(">>Per sequence quality scores"):
            print "[seed] All good, no seed needed."
            break
        #---------------------------------
    print '[parseFastqc] endstime %s.' % strftime(timeformat,\
            localtime())
#-------------End parseFastqc--------------------------------------

def estimateQuality(fastq):
    '''
    This function use fastqc to estimate the sequencing quality of
    fastq.
    '''
    print "[estimateQuality] Begins %s" \
        % strftime(timeformat, localtime())
    global totalReads
    if not os.path.exists(fastq):
        print >>sys.stder, "#No such file %s" & fastq
        sys.exit(1)
    #--------------------------------
    if totalReads == 0:
        totalReads = int(\
            os.popen('wc -l '+fastq).readline().split()[0])
    #------------------------------
    cmd = "fastqc -f fastq " + fastq
    try:
        os.system(cmd)
    except OSError, e:
        if e.errno == 32512:
            print >>sys.stderr, "No command fastqc found, please \
             check."
        else:
            print >>sys.stderr, "Unknown error code for %s" % cmd
        sys.exit(1)
    print "[%s] : Finish quality estamination. %s " % \
            (cmd, strftime(timeformat, localtime))
    if not seed:
        fastqc_data = os.path.join(\
            os.path.splitext(os.path.basename(fastq))[0] + '_fastqc',\
            'fastqc_data.txt')
        parseFastqc(fastqc_data)
#---------------------------------------------------------
def bwa(fastqc, index, parameter):
    '''
    This function does alignments.
    '''
    cmd = 'bwa aln -t 6 '
def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime   = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()

