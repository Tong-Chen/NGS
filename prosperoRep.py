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
from ctIO import readFasta, outputRep
import re

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) < 4:
        print >>sys.stderr, 'Using python %s pep prospero \
outputfile [overlap percentage]' % sys.argv[0]
        sys.exit(0)
    #---ori--------------------------------------------
    pat = re.compile(">.+?from (\d+) to (\d+).+?from (\d+) to (\d+) ")
    seqDict = readFasta(sys.argv[1])
    repDict = {}
    for line in open(sys.argv[2]):
        if line.startswith('using sequence1'):
            locus = line.strip().split()[-1]
            seq = seqDict[locus]
            repDict[locus] = []
        elif line[0] == '>':
            match = pat.match(line)
            tmpDict = {}
            pos1 = int(match.group(1))
            pos2 = int(match.group(2))
            pos3 = int(match.group(3))
            pos4 = int(match.group(4))
            if len(sys.argv) == 5:
                if (pos2-pos3+1.0)/(pos4-pos1+1.0) > \
                        float(sys.argv[4]):
                    continue
            tmpDict[(pos1, pos2)] = seq[pos1-1:pos2]  
            tmpDict[(pos3, pos4)] = seq[pos3-1:pos4]  
            repDict[locus].append(tmpDict)
    #------------------------------------------------
    outputRep(repDict, sys.argv[3])
if __name__ == '__main__':
    main()

