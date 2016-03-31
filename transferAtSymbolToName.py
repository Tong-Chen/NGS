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
import ctTEST
import re

pattern = re.compile(r'AT.G[0-9]{5}')

def readTransfer(file, flag = 'symbol', head = 1):
    if flag == 'symbol':
        symbolDict = {}
        flag = 1
    elif flag == 'locus':
        locusDict = {}
        flag = 2
    else:
        print >>sys.stderr, "Wrong flag, should be 'symbol' or 'locus'"
        sys.exit(1)
    #---------------------------------------------
    for line in open(file):
        if head:
            head -= 1
            continue
        #-------------------------
        line = line.rstrip()
        (locus, symbol, name) = line.split('\t')
        name = name.upper()
        if flag == 1:
            symbol = symbol.upper()
            if pattern.match(symbol):
                #print >>sys.stderr, symbol, name
                #assert symbol == name
                if symbol != name:
                    print >>sys.stderr, symbol, name
                    name = symbol
        if flag == 1:
            #There are duplicate symbols, before use list to save
            #names, later change to set.
            if symbol not in symbolDict:
                symbolDict[symbol] = set()
            symbolDict[symbol].add(name)
        elif flag == 2:
            #Do not know if there are duplicate locus
            if locus not in locusDict:
                locusDict[locus] = set()
            locusDict[locus].add(name)
        #-------------------------------------------
    #------------End reading------------------
    if flag == 1:
        if 0:
            ctTEST.ct_rdict(symbolDict)
            sys.exit()
        return symbolDict
    elif flag == 2:
        return locusDict
#-------------------------End readTransfer------------


def main():
    print >>sys.stderr, "This is used to translate gene locus or \
symbol of Arabidopsis to gene names. If symbol is not the same\
with name, let name equals symbol. If the symbol not in the list\
and match the pattern, print directly."
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #------------------------------------
    
    oriList = [line.strip().upper() for line in open(sys.argv[1])]

    transfer =\
    "/home/CT/server/project/project1/plant/Arabidopsis/TAIR.ID.SYMBOL.NAME.Translator"
    
    aDict = readTransfer(transfer)

    for locus in oriList:
        if locus in aDict:
            for item in aDict[locus]:
                if pattern.match(item):
                    print item
                else:
                    print >>sys.stderr, "-_-", locus, item
        else:
            if pattern.match(locus):
                print locus
            else:
                print >>sys.stderr, "-_-", locus
if __name__ == '__main__':
    main()

