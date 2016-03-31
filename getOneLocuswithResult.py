#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
import sys

def main():
    print >>sys.stderr, "This used to get the repetitions with setted\
locus. For example, to get the real protein repetition\nPrint the\
 result to the screen"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s locusfile resultfile' % sys.argv[0]
        sys.exit(0)
    
    rDict = {}
    for line in open(sys.argv[2]):
        if line[0] == '>':
            locus = line[1:-1]
            rDict[locus] = []
        else:
            rDict[locus].append(line)
    #----------------------------------------
    keya = rDict.keys()
    for locus in open(sys.argv[1]):
        locus = locus[:-1]
        for i in keya:
            if i.find(locus) != -1:
                print '>%s' % i
                for i33 in rDict[i]:
                    print i33,
    #------------------------------------------
if __name__ == '__main__':
    main()

