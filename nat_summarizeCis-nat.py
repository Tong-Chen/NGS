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
    This is designed to summarize the statistics of cis-Natural
    antisense transcripts.

Input file:

for -i:

    (   
        Normally input file is generated using command 
        #intersectBed -a tr.bed -b tr.bed -wo -S >output
    )

#The name of each bed region has special formats as explained below:
    f.I.1: 
        f-represents gene name (no '@' in gene name); 
        I-represents special parts of gene, here is intron.
        C-represents special parts of gene, here is coding-exon.
        1-represents the first intron or coding-exon.
        
        One can use other different alphabets or words (no dot '.' allowed) to
        represent intron or exon.

chr1    9550913 9621173 b       1       -       chr1    9602747 9603392 f@I@1   0       +       645
chr1    9550913 9621173 b       1       -       chr1    9603392 9603555 f@C@2   0       +       163
chr1    9550913 9621173 b       1       -       chr1    9603392 9603555 g@C@3   0       +       163
chr1    9550913 9621173 b       1       -       chr1    9603555 9609916 f@I@2   0       +       6361
chr1    9550913 9621173 b       1       -       chr1    9603555	9609916 g@I@3   0       +       6361
chr1    9550913 9621173 b       1       -       chr1    9609916	9609978 f@C@3   0       +       62
chr1    9550913 9621173 b       1       -       chr1    9609916	9609978 g@C@4   0       +       62

############################# un used #####################
###for -g:

###The name for each bed region in <transcript-position-file> should be
###consistent with the transcript name in file given to <-i>.

###Normally this file is part of input file for intersectBed.
############################# un used #####################
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

verbose=0

def fprint(content):
    print >>sys.stderr, json_dumps(content,indent=1)

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
        metavar="FILEIN", help="Normally the output of <intersectBed -a tr.bed -b tr.bed -wo -S >output>")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def get_overlap(posL):
    '''
    posL = [left_start, left_end, right_start, right_end]
    '''
    left = posL[0] if posL[0]>posL[2] else posL[2]
    right = posL[1] if posL[1]<posL[3] else posL[3]
    if left <= right:
        return right - left
    else:
        return -1
#--------------------------------------------------------------------

def determineMaxPercent(valueL):
    '''
    valueL = [left_start, left_end, 
        right_start, right_end, left_strand, right_strand, 
        overlap, match_type, chr]
    '''
    left_start = valueL[0]
    left_end   = valueL[1]
    left_len = left_end - left_start
    right_start = valueL[2]
    right_end   = valueL[3]
    right_len = right_end - right_start
    minLen = left_len if left_len < right_len else right_len
    return float(valueL[6]) / minLen
    
#---------------------END determineOverlap_pos-----------------
#-----overlap-------------------------------
def determineOverlap_pos(valueL):
    '''
    valueL = [left_start, left_end, 
        right_start, right_end, left_strand, right_strand, 
        overlap, match_type, chr]
    '''
    ####Only for test-----------------------------
    if verbose:
        overlap = get_overlap(valueL[0:4])
        assert overlap == valueL[6]
    #-------determine head-head tail-tail contain-----
    left_start = valueL[0]
    left_end   = valueL[1]
    right_start = valueL[2]
    right_end   = valueL[3]
    left_strand = valueL[4]

    if right_end > left_end > right_start > left_start:
        if left_strand == '+':
            type = 'tail-2-tail'
        else:
            type = 'head-2-head'
    elif left_end > right_end > left_start > right_start:
        if left_strand == '+':
            type = 'head-2-head'
        else:
            type = 'tail-2-tail'
    else:
        type = 'full-in'
    #-------------------------------------------
    return type
#-----------------------------------------------
def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    #trpos = options.trpos
    global verbose
    verbose = options.verbose
    debug = options.debug
    symbol = '@'
    #------trPos dict-------------------
    #trposD = {}
    #for line in open(trpos):
    #    lineL = line.split()
    #    tr = lineL[3]
    #    if tr in trposD:
    #        print >>sts.stderr, "Duplicate tr %s" % tr
    #    #-------------------------------------------
    #    trposD[tr] = [int(lineL[1]), int(lineL[2]), lineL[5]] 
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    for line in fh:
        lineL = line.split('\t')
        chr = lineL[0]
        left_name = lineL[3]
        left_nameL = left_name.split(symbol)
        len_left_nameL = len(left_nameL)
        left_type = 'tr'
        left_tr = left_nameL[0]
        if len_left_nameL >= 2:
            left_type = left_nameL[1]
        left_start = int(lineL[1])
        left_end = int(lineL[2])
        left_strand = lineL[5]
        #if lineL[4] != '0':
        #    trposD[left_name] = [left_start, left_end, left_strand]
        chr = lineL[0]
        right_name = lineL[9]
        right_nameL = right_name.split(symbol)
        len_right_nameL = len(right_nameL)
        right_type = 'tr'
        right_tr = right_nameL[0]
        if len_right_nameL >= 2:
            right_type = right_nameL[1]
        right_start = int(lineL[7])
        right_end = int(lineL[8])
        right_strand = lineL[11]
        #if lineL[10] != '0':
        #    trposD[right_name] = [right_start, right_end, right_strand]
        overlap = int(lineL[12])
        key = (left_name, right_name)
        match_type = [[left_type], [right_type]]
        if key not in aDict:
            aDict[key] = [left_start, left_end, 
                right_start, right_end, left_strand, right_strand, 
                overlap, match_type, chr]
        else:
            if left_start < aDict[key][0]:
                aDict[key][0] = left_start
            if left_end > aDict[key][1]:
                aDict[key][1] = left_end
            if right_start < aDict[key][2]:
                aDict[key][2] = right_start
            if right_end > aDict[key][3]:
                aDict[key][3] = right_end
            #--------------------------------------
            aDict[key][6] = get_overlap(aDict[key][0:4]) 
            aDict[key][7][0].append(left_type)
            aDict[key][7][1].append(right_type)
        #-------------------------------------------------
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    if verbose == '1':
        for key, valueL in aDict.items():
            print >>sys.stderr, key
            print >>sys.stderr, valueL
    ###-----------------------------------------
    duplicateD = {}
    '''
    valueL = [left_start, left_end, right_start, right_end, 
              left_strand, right_strand, overlap, 
              match_type, chr]
    '''
    header = ['left_name', 'left_pos', 'right_pos', 'left_strand',
            'right_name', 'right_pos', 'right_pos', 'right_strand', 
            'overlap', 'overlap_percent', 'match_type', 'chr',
            'pair_type']
    print '\t'.join(header)
    for keyT, valueL in aDict.items():
        if keyT in duplicateD:
            continue
        pair_pos = determineOverlap_pos(valueL)
        maxPercent = determineMaxPercent(valueL)
        match_type = '@'.join(['-'.join(sorted(pair)) for pair in
            sorted(valueL[7])])
        print "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%.2f\t%s\t%s\t%s" % \
                (keyT[0], valueL[0], valueL[1], valueL[4],
                keyT[1], valueL[2], valueL[3], valueL[5], valueL[6],
                maxPercent, match_type, valueL[8], pair_pos)
        duplicateD[(keyT[1], keyT[0])] = 1
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


