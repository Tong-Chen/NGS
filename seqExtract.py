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

def readGenome(genome):
    aDict = {}
    key = ''
    for line in open(genome):
        if line[0] == '>':
            if key:
                aDict[key] = ''.join(alist)
            #--------------------------------
            key = line[1:-1]
            alist = []
        else:
            alist.append(line.upper().strip())
    #-------------------------------------------
    aDict[key] = ''.join(alist)
    return aDict
#--------------------------------

def getSeq(bed,  aDict):
    '''
    bed: 0-based,  half-open
    '''
    output = os.path.split(bed)[1] + '.fa'
    fh = open(output, 'w')
    for line in open(bed):
        lineL = line.strip().split('\t')
        chr, start, end, name = lineL[:4] 
        start = int(start)
        end   = int(end)
        print >>fh, '>%s\n%s' % (name,  aDict[chr][start:end])
    fh.close()
#-------------------------------------
#------------------------------------------
def main():
    if len(sys.argv) < 3:
        print >>sys.stderr, "Get sequence of a bed file from genome."
        print >>sys.stderr, "Print the result to file." 
        print >>sys.stderr, 'Using python %s genome \
pos_bed1[,pos_bed2...](at least four cols, the 4th col is the name of\
 fasta sequence)' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------------------
    genome = sys.argv[1]
    bedL = sys.argv[2:]
    aDict = readGenome(genome)
    for bed in bedL:
        getSeq(bed,  aDict)

if __name__ == '__main__':
    main()

