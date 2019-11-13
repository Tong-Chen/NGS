#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, with_statement
'''
Copyright 2013, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Functional description:
    This is used to summarize the length of FASTQ sequences.
'''

import sys
import gzip
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
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--fastq-name", dest="fastq",
        metavar="FILEIN", help="FASTQ file name for isequencing \
reads without fq.gz like T1.")
    parser.add_option("-t", "--fastq-type", dest="fastq_type",
        metavar="PE/SE", help="FASTQ file type <PE (pair-end default) or <SE>.")
    parser.add_option("-s", "--suffix", dest="suffix",
        default=".fq.gz", metavar="SUFFIX", help="Normally <.fq.gz (default)> or <.fq>.")
    parser.add_option("-a", "--auto-dertermine-len", dest="auto_len",
        default=1, type='int', help="Given <1> to let the program select \
best length. Default <1> means determining length parameter.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.fastq1 != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def determineMostLen(aDict):
    '''
    aDict = {(30, 30):4, (30, 29):1}
    '''
    len_pairL = aDict.keys()
    for len_left, len_right in len_pairL:
        if len_left > len_right > 0:
            count = aDict.pop((len_left, len_right))
            len_left = len_right
            len_pair = (len_left,len_right)
            aDict[len_pair] = aDict.get(len_pair, 0)+count
        elif 0< len_left < len_right:
            count = aDict.pop((len_left, len_right))
            len_right = len_left
            len_pair = (len_left,len_right)
            aDict[len_pair] = aDict.get(len_pair, 0)+count
    #----------------------------------------
    pairL = aDict.keys()
    pairL.sort(reverse=True)
    len_d = len(pairL)
    count = 0
    maxLen = pairL[0][0]
    allBase = 2 * sum([key[0]*value for key,value in aDict.items()])
    baseD = {}
    for i in range(len_d):
        currentP = pairL[i]
        count += aDict[currentP]
        currentL = currentP[0]
        savePercent = currentL*2*count / allBase
        baseD[currentL] = [currentL, savePercent, savePercent*currentL]
    #--------------------------------
    valueL = baseD.values()
    valueL.sort(key=lambda x:x[1], reverse=True)
    for value in valueL:
        print >>sys.stderr,value
    return valueL[0][0]
#-----------------------------------

def getLen(fq1, fq2, fastq1, fastq2):
    aDict = {}
    lenD = {}
    i = 0
    for line in fq1:
        i += 1
        line2 = fq2.readline() if fq2 else ''
        if i % 4 == 2:
            len_1 = len(line.strip())
            len_2 = len(line2.strip())
            len_pair = (len_1, len_2)
            aDict[len_pair] = aDict.get(len_pair, 0)+1
    #----------END fq_h--------------------------
    print >>sys.stderr, "Length distribution for %s %s" % (fastq1, fastq2)
    len_iK = aDict.keys()
    len_iK.sort(key=lambda x: sum(x), reverse=True)
    for len_pair in len_iK:
        print >>sys.stderr, "(%d, %d)\t%d" % (len_pair[0], len_pair[1], aDict[len_pair])
    #---------------------------------------------
    return determineMostLen(aDict)
#---------------------------------------
    
#---------------------------

#---------------------------------------
def trim(fq1, fq2, length, fq1_out, fq2_out):
    count = 0
    while 1:
        count += 1
        name1  = fq1.readline()
        if not name1:
            break
        seq1   = fq1.readline()
        len1 = len(seq1)
        third1 = fq1.readline()
        quan1  = fq1.readline()
        shorter = len1
        if fq2:
            name2  = fq2.readline()
            seq2   = fq2.readline()
            len2 = len(seq2)
            third2 = fq2.readline()
            quan2  = fq2.readline()
            shorter = len1 if len1 < len2 else len2
        if shorter < length:
            continue
        #----------------------------
        start = 0
        sec_count = 0
        for i in range(length, shorter, length):
            sec_count += 1
            name1 = '@%s_%s' % (str(count), str(sec_count)) 
            newseq1 = seq1[start:i]
            third1 = '+'
            newquan1 = quan1[start:i]
            print >>fq1_out, '\n'.join([name1, newseq1, third1, newquan1])
            if fq2:
                name2 = '@%s_%s' % (str(count), str(sec_count)) 
                newseq2 = seq2[start:i]
                third2 = '+'
                newquan2 = quan2[start:i]
                print >>fq2_out, '\n'.join([name2, newseq2, third2, newquan2])
            #--------------------------------------
            start = i
        #------------------------------------------
#----------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    suffix = options.suffix
    fastq1 = options.fastq1 + suffix
    if options.fastq2:
        fastq2 = options.fastq2 + suffix
    else:
        fastq2 = ''
    fq1 = fq2 = ''
    auto_len = options.auto_len
    gz = 0
    if fastq1.endswith('.gz'):
        gz = 1
    #---------------------------------
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if gz:
        fq1 = gzip.open(fastq1, 'rb')
        fq2 = gzip.open(fastq2, 'rb') if fastq2 else ''
    else:
        fq1 = open(fastq1, 'r')
        fq2 = open(fastq2, 'r') if fastq2 else ''
    #-----------------------------------
    if auto_len:
        length = getLen(fq1, fq2, fastq1, fastq2)
        print length
    fq1.close()
    if fastq2:
        fq2.close()
    #----------------------------------------------------
#    if gz:
#        fq1 = gzip.open(fastq1, 'rb')
#        if length:
#            fq1_out_file = options.fastq1+'.trim'+str(length) + suffix
#            fq1_out = gzip.open(fq1_out_file, 'wb')
#        if fastq2:
#            fq2 = gzip.open(fastq2, 'rb')
#            if length:
#                fq2_out_file = options.fastq2+'.trim'+str(length) + suffix
#                fq2_out = gzip.open(fq2_out_file, 'wb')
#        #------------------------------------------
#    else:
#        fq1 = open(fastq1, 'r')
#        if length:
#            fq1_out_file = options.fastq1+'.trim'+str(length) + suffix
#            fq1_out = open(fq1_out_file, 'w')
#        if fastq2:
#            fq2 = open(fastq2, 'r')
#            if length:
#                fq2_out_file = options.fastq2+'.trim'+str(length) + suffix
#                fq2_out = open(fq2_out_file, 'w')
#        #------------------------------------------
#    #-----------------------------------
#    if fastq2:
#        trim(fq1, fq2, length, fq1_out, fq2_out)
#    else:
#        fq2 = fq2_out = ""
#        trim(fq1, fq2, length, fq1_out, fq2_out)
#    #----------------
#    fq1.close()
#    fq1_out.close()
#    if fastq2:
#        fq2.close()
#        fq2_out.close()
#    ###--------multi-process------------------
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


