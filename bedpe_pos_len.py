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
    This is designed to extract the start and end position of PE reads (using the start of first read and end of second read), length, count.
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
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="BEDPE file")
    parser.add_option("-s", "--suppress-short-regions", dest="short",
        default=100, type='int', 
        help="Suppress regions shorted than given value. This is mainly designed to exclude soft clipped reads.")
    parser.add_option("-g", "--genome", dest="genome",
        help="Genome sequence in fasta format.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def readGenome(genome):
    genomeD = {}
    key = ''
    seqL = ''
    for line in open(genome):
        if line[0] == '>':
            if seqL and key:
                genomeD[key] = ''.join(seqL)
            key = line[1:-1]
            assert key not in genomeD, key
            seqL = []
        else:
            seqL.append(line.strip())
    #--------------------------------------
    if seqL and key:
        genomeD[key] = ''.join(seqL)
    return genomeD
#--------------------------------
def getSeq(genomeD, key):
    chr, start, end, length = key.split('\t')
    start = int(start)
    end = int(end)
    return genomeD[chr][start:end].upper()
#--------------------------------
#def actgPercent(seq):
#    A = seq.count('A')
#    T = seq.count('T')


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    verbose = options.verbose
    short = options.short
    genome = options.genome
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    for line in fh:
        lineL = line.strip().split('\t')    
        assert lineL[0] == lineL[3], lineL 
        #assert lineL[1] <= lineL[4], lineL 
        #assert lineL[2] <= lineL[5], lineL 
        assert lineL[7] != lineL[8], lineL
        posL = [int(lineL[1]), int(lineL[2]), int(lineL[4]), int(lineL[5])]
        posL.sort()
        start = posL[0]
        end   = posL[-1]
        assert start < end, lineL
        length = end-start
        if length < short:
            continue
        read = '{}\t{}\t{}\t{}'.format(lineL[0], start, end, length)
        aDict[read] = aDict.get(read, 0) + 1
        
    #-------------END reading file----------
    genomeD = readGenome(genome)

    total_count = sum(aDict.values()) / (10**6)
    resultD = {}
    print "chr\tstart_0_based\tend_0_based_not_include\tLength\tCount\tSeq\tCounts_per_million\tA_percent\tT_percent\tC_percent\tG_percent\tCpG_percent\tCpG_observe_to_expect_ratio"
    keyL = aDict.keys()
    keyL.sort(key=lambda x: aDict[x], reverse=True)
    for key in keyL:
        seq = getSeq(genomeD, key)
        len_seq = len(seq)
        A = "%.2f" % (100 * seq.count('A') / len_seq)
        T = "%.2f" % (100 * seq.count('T') / len_seq)
        C = seq.count('C')
        G = seq.count('G')
        CG = seq.count('CG')
        '''
        The CpG count is the number of CG dinucleotides in the island. The Percentage CpG is the ratio of CpG nucleotide bases (twice the CpG count) to the length. The ratio of observed to expected CpG is calculated according to the formula (cited in Gardiner-Garden et al. (1987)):
        
                Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G)
                
                where N = length of sequence.
        '''
        CG_per = "%.2f" % (100 * CG*2/len_seq)
        if C*G == 0:
            CG_obvsexpect = 0
        else:
            CG_obvsexpect = "%.2f" % (100 * CG / (C*G) * len_seq)
        C = "%.2f" % (100 * C / len_seq)
        G = "%.2f" % (100 * G / len_seq)
        print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format\
                (key, aDict[key], seq, aDict[key]/total_count, A, T, C, G, CG_per, CG_obvsexpect) 
    #----close file handle for files-----
    if file != '-':
        fh.close()
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


