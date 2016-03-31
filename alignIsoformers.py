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
import ctIO
import os

def main():
    print >>sys.stderr, "Print the result to current directory"
    lenS = len(sys.argv)
    if lenS != 3 and lenS != 2:
        print >>sys.stderr, 'Using python %s seqfile(for.c) [locus](optional)' % sys.argv[0]
        sys.exit(0)
    #--------------------------------
    seqDict = {}
    ctIO.readseq(sys.argv[1], seqDict)
    locusDict = {}
    if lenS == 2:
        locusL = seqDict.keys()
        locusL.sort()
    elif lenS == 3:
        locusL = [line.strip() for line in open(sys.argv[2])]
    for line in locusL:
        line = line.strip()
        dot = line.find('.')
        locus = line[:dot]
        after = line[dot:]
        if locus in locusDict:
            locusDict[locus].append(after)
        else:
            locusDict[locus] = [after]
    #--------------------------------
    fhl = open('locusHasIsoforms', 'w')
    for locus, afterL in locusDict.items():
        if len(afterL) < 2:
            continue
        filename = locus+'.fa'
        
        fh = open(filename, 'w')
        for after in afterL:
            newLoc = locus+after
            print >>fhl, newLoc
            print >>fh, '>%s' % newLoc
            print >>fh, seqDict[newLoc]
        fh.close()
        cmd = 't_coffee ' + filename
        os.system(cmd)
    #---------------------------
    fhl.close()
if __name__ == '__main__':
    main()

