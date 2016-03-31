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
    parser.add_option("-i", "--fastq1", dest="fastq1",
        metavar="FILEIN", help="FASTQ file name without \
fq.gz like T1_1.")
    parser.add_option("-j", "--fastq2", dest="fastq2",
        metavar="FILEIN", help="FASTQ file name without \
fq.gz like T1_2.")
    parser.add_option("-s", "--suffix", dest="suffix",
        metavar="SUFFIX", help="Normally .fq.gz or .fq.")
    parser.add_option("-l", "--length", dest="length",
        default=0, help="Reads smaller than this length will be removed. \
Longer than this length will be trimmed. Multiple this length \
will be split. Default 0 means no trim needed. ")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.fastq1 != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def getLen(fq_h, fastq):
    aDict = {}
    lenD = {}
    i = 0
    for line in fq_h:
        i += 1
        if i % 4 == 2:
            len_i = len(line.strip())
            if len_i not in aDict:
                aDict[len_i] = len_i
            else:
                aDict[len_i] += len_i
    #----------END fq_h--------------------------
    print "Length distribution for %s" % fastq
    len_iK = aDict.keys()
    len_iK.sort()
    for len_i in len_iK:
        print "%d\t%d" % (len_i, aDict[len_i])
    #---------------------------------------------
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
    fastq2 = options.fastq2 + suffix
    fq1 = fq2 = ''
    length = int(options.length)
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
        if length:
            fq1_out_file = options.fastq1+'.trim'+str(length) + suffix
            fq1_out = gzip.open(fq1_out_file, 'wb')
        if fastq2:
            fq2 = gzip.open(fastq2, 'rb')
            if length:
                fq2_out_file = options.fastq2+'.trim'+str(length) + suffix
                fq2_out = gzip.open(fq2_out_file, 'wb')
        #------------------------------------------
    else:
        fq1 = open(fastq1, 'r')
        if length:
            fq1_out_file = options.fastq1+'.trim'+str(length) + suffix
            fq1_out = open(fq1_out_file, 'w')
        if fastq2:
            fq2 = open(fastq2, 'r')
            if length:
                fq2_out_file = options.fastq2+'.trim'+str(length) + suffix
                fq2_out = open(fq2_out_file, 'w')
        #------------------------------------------
    #-----------------------------------
    if not length:
        getLen(fq1, fastq1)
        if fastq2:
            getLen(fq2, fastq2)
    else:
        if fastq2:
            trim(fq1, fq2, length, fq1_out, fq2_out)
        else:
            fq2 = fq2_out = ""
            trim(fq1, fq2, length, fq1_out, fq2_out)
    #----------------
    fq1.close()
    fq1_out.close()
    if fastq2:
        fq2.close()
        fq2_out.close()
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


