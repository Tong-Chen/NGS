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

def purify(file):
    aDict = {}
    head = 1
    for line in open(file):
        if head:
            keyL = line.split('\t')
            lenKeyL = len(keyL)
            for item in keyL:
                aDict[item] = []
            head -= 1
        else:
            valueL = line.split('\t')
            for i in range(lenKeyL):
                aDict[keyL[i]].append(valueL[i])
    #-------------------------------------------
    newKeyL = [keyL[0]]
    for key in keyL[1:]:
        newvalue = [int(item) for item in aDict[key]]
        if sum(newvalue) == 0:
            del aDict[key]
        else:
            newKeyL.append(key)
    #--------------------------------------------
    newfile = file+'.Purify'
    lennewKeyL = len(newKeyL)
    lenvalueL = len(aDict[newKeyL[0]])
    fh = open(newfile, 'w')
    print >>fh, '\t'.join(newKeyL)
    for j in range(lenvalueL):
        lineL = []
        for i in range(lennewKeyL):
            lineL.append((aDict[newKeyL[i]])[j])
        #-----------------
        print >>fh, '\t'.join(lineL)

    #---------------------------------------------
    fh.close()
#-------------------------------------------
def main():
    print >>sys.stderr, "Purify the result, delete the column full\
with 0. And output the result to prefix.Purify"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s prefix' % sys.argv[0]
        sys.exit(0)
    #----------------------------------
    prefix = sys.argv[1]
    suffix = ['GoBpSymbol','GoCcSymbol','GoMfSymbol','InterproSymbol',\
        'PfamSymbol', 'GoBpDescrip','GoCcDescrip','GoMfDescrip',\
        'InterproDescrip','PfamDescrip']
    for i in suffix:
        purify(prefix+'.'+i)
if __name__ == '__main__':
    main()

