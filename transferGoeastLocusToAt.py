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

def readTransfer(transfer): 
    head = 1
    aDict = {}
    for line in open(transfer):
        if head:
            head -= 1
            continue
        #-----------------------------
        locus, symbol, name = line.strip().split('\t')
        if locus not in aDict:
            aDict[locus] = ''
        else:
            print >>sys.stderr, "Duplicated locus %s" % locus
        if symbol != name:
            if symbol.startswith('AT') and len(symbol) == 9:
                aDict[locus] = symbol
            elif name.startswith('AT') and len(name) == 9:
                aDict[locus] = name
            else:
                print >>sys.stderr, '%s\t%s\t%s' % (locus, symbol,
                    name)
            #------------------------------
        else:
            aDict[locus] = name
    return aDict
#----------------------------------------------

def main():
    print >>sys.stderr, "This is used to transfer the arabidopsis locus\
in GOEAST result to AT started name. Print to standard output."
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s goeastOutput' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------
    transfer = \
            "/home/CT/server/project/project1/plant/Arabidopsis/TAIR.ID.SYMBOL.NAME.Translator"
    aDict = readTransfer(transfer)
    head = 1
    for line in open(sys.argv[1]):
        if head:
            print line,
            head -= 1
            continue
        #------------------------------
        lineL = line.strip().rsplit('\t')
        front = '\t'.join(lineL[0:7])
        middle = lineL[7]
        end = '\t'.join(lineL[8:])
        middleL = middle.split(' // ')
        middleLAt = [aDict[locus] for locus in middleL]
        print '%s\t%s\t%s' % (front, ' // '.join(middleLAt), end)
if __name__ == '__main__':
    main()
