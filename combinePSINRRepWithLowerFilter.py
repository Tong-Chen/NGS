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
from ctIO import readFasta

def readSubjS(subjS):
    '''
    '''
    subjSDict = {}
    for line in open(subjS):
        if line[0] == '>':
            assert(line[:3] == ">gi")
            gi = line.split('|')[1]
            if gi not in subjSDict:
                subjSDict[gi] = ''
            else:
                print >>sys.stderr, 'Duplicate gi', gi, 'in', subjS
        else:
            subjSDict[gi] += line.strip()
        #------End one line--------------------
    #----------End one file-------------------
    return subjSDict
#-------------------------------------------------

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 4:
        print >>sys.stderr, 'Using python %s filename subjS atseq' % sys.argv[0]
        sys.exit(0)
    #---------------------------------------------------
    subjSDict = readSubjS(sys.argv[2])
    atDict = readFasta(sys.argv[3])
    at = 1
    for line in open(sys.argv[1]):
        if line[0] == '=':
            group = line[1:].split()[1]
            at = 1  #label the following locus is Arabidopsis
        elif line[0] == '>':
            if at:
                locus = (line[1:].rsplit('.', 1))[0]
                seq = atDict[locus]
                at = 0
            else:
                locus = line[1:-1]
                seq = subjSDict[locus]
            #--------------------------------
        else:

if __name__ == '__main__':
    main()

