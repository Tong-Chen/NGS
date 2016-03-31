#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
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
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #--------------------------------------------------------------
    print '>locus\tanyMutPart\tsecMutPart\tanyMutAll\tsecMutAll'
    for line in open(sys.argv[1]):
        if line[0] == '>':
            print line,
        else:
            tmpdict = {}
            line = line.rstrip()
            tmplist = line.rstrip('#').split('#')
            seqlist = [x.split(':')[0] for x in tmplist]
            repNum = len(seqlist)
            if repNum > 4:
                repNum = 4
            del tmplist
            length = len(seqlist[0])
            for seq in seqlist:
                tmplen = len(seq)
                if length > tmplen:
                    length = tmplen
                    print 'A',
            for pos in range(length):
                tmpdict[pos] = set()
                for seq in seqlist:
                    tmpdict[pos].add(seq[pos])
                #----------------------------------------------
            #-------------------------------------------------
            anyMutPart = 0 #any position have a mutant
            secMutPart = 0 # second position have a mutant
            anyMutAll = 0 #any position all mutant
            secMutAll = 0 #second position all mutant
            i = 0
            for value in tmpdict.values():
                count = len(value)
                i += 1
                if count > 1:
                    anyMutPart += 1
                    if i % 3 == 2:
                        secMutPart += 1
                if count == repNum:
                    anyMutAll += 1
                    if i % 3 == 2:
                        secMutAll += 1
                #-------------------------------
            del tmpdict
            #---------------------------------------------
            print '\t%.1f\t%.1f\t%.1f\t%.1f' % \
                (anyMutPart/length, secMutPart/length,\
                anyMutAll/length, secMutAll/length)
                


if __name__ == '__main__':
    main()

