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

Four columns in output file (from left to right) represent the context
for A, C, G, U ordinal.

Score file

15  3.87880514112849
16  4.66726250149269
17  5.45571986185696
18  6.24417722222123
19  7.03263458258551
20  7.82109194294978
21  8.60954930331405
22  9.39800666367831
23  10.1864640240426
24  10.9749213844069
25  11.7633787447711
26  12.5518361051354
27  13.3402934654997
28  14.1287508258640
29  14.9172081862283

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
        metavar="FILEIN", help="miRNA mature sequences, FASTA format\
, all sequences in one line.")
    parser.add_option("-s", "--score-file", dest="score",
        metavar="score-file", default="/home/chentong/home/soft/homer/score.txt", 
        help="tab splitted lines with the first \
column as length and the second column as score. \
There is a default file for this parameter.")
    parser.add_option("-n", "--value-for-nonSeed", dest="non_seed",
        metavar="0.55", default=0.55, help="A number to represent the weight of \
given nucleotide.")
    parser.add_option("-t", "--value-for-seed", dest="seed",
        metavar="0.997", default=0.997, help="A number to represent the weight of \
given nucleotide.")
    parser.add_option("-c", "--complement", dest="complement",
        default=1, help="Do you want to get the complementary sequence \
for miRNAs, default TRUE. Accept '0' to turn it off. ")
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
    score = options.score
    nonSeedValue = float(options.non_seed)
    seedValue = float(options.seed)
    complement = int(options.complement)
    verbose = options.verbose
    debug = options.debug
    dict = {'A':'U','G':'C','U':'A','C':'G','a':'u','g':'c','u':'a','c':'g'}
    ntPos = {'A':0, 'C':1, 'G':2, 'U':3}
    ntL = ['A', 'C', 'G', 'U']
    #-----------------------------------
    scoreD = {}
    for line in open(score):
        key, value = line.split()
        scoreD[key] = value
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    for line in fh:
        if line[0] == '>':
            key = line.split()[0][1:]
        else:
            matrixL = []
            seq = line.strip()
            lenseq = len(seq)
            score = scoreD[str(lenseq)]
            #print complement
            if complement != 0:
                #print 'here'
                comseq = [dict[i] for i in line.strip()]
            else:
                comseq = [i for i in line.strip()]
            #sys.exit(2)
            #----------------------
            cnt = 0
            for nt in comseq:
                tmpD = {'A':0, 'C':0, 'G':0, 'U':0}
                cnt += 1
                if 2 <= cnt <= 8: 
                    tmpD[nt] = seedValue
                else:
                    tmpD[nt] = nonSeedValue
                #-------------------------
                other = (1-tmpD[nt]) / 3
                for key89,v89 in tmpD.items():
                    if v89 == 0:
                        tmpD[key89] = other
                #------------------------------
                matrixL.append('\t'.join([str(tmpD[i93]) for i93 in
                    ntL]))
            #------------------------
            matrixL.reverse()
            print ">%s\t%s\t%s" % (seq,key,score)
            print "\n".join(matrixL)
        #---------------------------------------
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
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



