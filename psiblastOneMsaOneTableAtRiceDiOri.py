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

def readSimple(file):
    '''
    aDict = {'AT':[LOC1, LOC2]}
    '''
    aDict = {}
    for line in open(file):
        if line[0] == '>' and line[-3:-1] != ' 0':
            locus = line[1:].split()[0]
            aDict[locus] = []
        elif line[0] != '>':
            aDict[locus] = line.strip().split()
    #----------------------------------
    return aDict
#-----------------------
def compare(atDict, sativaDict):
    for key, valueL in atDict.items():
        for item in valueL:
            if item in sativaDict:
                if key in sativaDict[item]:
                    print "%s\t%s" % (key, item)

def main():
    print >>sys.stderr, "Print the result to screen"
    print >>sys.stderr, "This used to parse the result produced by \
psiblastOneMsaOneTableAtRice.py. Specially this one is used to detect \
atrice comparision analysis and get di-orientation results. "
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s atsimple sativasimple' % sys.argv[0]
        sys.exit(0)
    #--------------------------------------
    atDict = readSimple(sys.argv[1])
    sativaDict = readSimple(sys.argv[2])
    compare(atDict, sativaDict)

if __name__ == '__main__':
    main()

