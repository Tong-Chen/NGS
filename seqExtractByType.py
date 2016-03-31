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
from ctFunc import find_all

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

def getSeq(typeL, aDict):
    '''
    bed: 0-based,  half-open
    '''
    #output = '.'.join([prefix, '.'.join(typeL), '.bed'])
    for chr, seq in aDict.items():
        #print seq
        for alpha in typeL:
            for i in find_all(seq, alpha):
                print "%s\t%d\t%d\t%s" % (chr, i, i+1, alpha)
    #------------------------------------------------ 
#-------------------------------------

def main():
    print >>sys.stderr, "This outputs a bed file to give the position of inputted alphabets"
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) < 3:
        print >>sys.stderr, 'Using python %s genome \
alphabet1[,alphabet2...]' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------------------
    genome = sys.argv[1]
    typeL = sys.argv[2:]
    aDict = readGenome(genome)
    getSeq(typeL, aDict)
if __name__ == '__main__':
    main()

