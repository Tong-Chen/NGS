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
#from re import compile

def addDict(adict, key, value):
    if key not in adict:
        adict[key] = set()
    adict[key].add(value)

def output(aDict):
    for key, value in aDict.items():
        fh = open(key, 'w')
        print >>fh, '\n'.join(value)
        fh.close()

def main():
    print >>sys.stderr, "Transfer aracyc_pathways to seperated subfiles"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------------------
    #pat = compile('(.*?) +(\d\.\d\.\d\.[^ ]*) +(.*?) +([^ ]+$)')
    firstDict = {} #1,2,3,4,5,6
    secondDict = {} #1.1, 1.2, .....6.9
    head = 1
    for line in open(sys.argv[1]):
        if head:
            head -= 1
            continue
        #matchO = pat.match(line.strip())
        #print '*',line,
        #if not matchO:
        #    print '*',line,
        #    continue
        #print '\n'.join((matchO.group(1),matchO.group(2),matchO.group(3),\
        #    matchO.group(4)))
        #reaction = matchO.group(2)
        #gene = matchO.group(4)
        lineL = line.strip().split('\t')
        reaction = lineL[1]
        gene = lineL[3]
        if gene != 'unknown':
            if not gene.startswith('AT'):
                print gene
                continue
            #print lineL
            #print reaction
            #this used to detect how many reactions not numbered
            #try:
            #    assert reaction.count('.') == 3
            #except AssertionError:
            #    print lineL
            #    print reaction
            #    continue
            if reaction.count('.') != 3:
                print '*', line
                addDict(firstDict, '7', gene)
            else:
                reactionL = reaction.split('.',2)
                firstK = reactionL[0]
                secondK = '.'.join(reactionL[0:2])
                addDict(firstDict, firstK, gene)
                addDict(secondDict, secondK, gene)
        else:
            print '#', line
    #--------------------------------------------------
    output(firstDict)
    output(secondDict)
if __name__ == '__main__':
    main()

