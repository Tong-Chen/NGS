#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, with_statement
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
import sys

def rdict(conDict):
    for locus, valueD in conDict.items():
        print locus
        for key, value in valueD.items():
            print '==', key

def conservationRate(conDict, numHomolog):
    #to deal with no hit, usually this will not happen.
    if len(conDict) == 0:
        print >>sys.stderr, 'No hit'
        print 0
        return 0
    for locus, valueD in conDict.items():
        valueDK = valueD.keys()
        lenvalueDK = len(valueDK)
        if lenvalueDK == 1:
            if locus.find(valueDK[0]) != -1:
                print 0
            else:
                print >>sys.stderr, 'No self hit', locus
                print "%.2f" % (1 / numHomolog)
        #actually this will not happen
        elif lenvalueDK == 0:
            print >>sys.stderr, 'No hit', locus
            print 0
        else:
            hitself = 0
            for hit in valueDK:
                if locus.find(hit) != -1:
                    hitself = 1
                    break
            if hitself:
                print "%.2f" % ((lenvalueDK - 1) / numHomolog)
            else:
                print >>sys.stderr, 'No self hit', locus
                print "%.2f" % (lenvalueDK / numHomolog)
        #--------End of else, one query----------------
    #------------End all query----------------------------
#------------End function-------------------------------------

            

def main():
    if len(sys.argv) != 5:
        print >>sys.stderr, 'Using python %s table subject identity[70] \
        evalue[0.01]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------------------------
    '''
    # Fields: query id, subject id, % identity, alignment length,
    # mismatches, gap opens, q. start, q. end, s. start, s. end, evalue,
    # bit score
    AT3G22120.1.6.169:189:219:229   AT3G22120.1     100.00  10      0
    0       1       10      169     178     7e-07   35.0

    conDict = {locus: {hitlocus:([],)}}
    '''
    #patch for select identity and pvalue. 20110930
    #identity = 0.7 pvalue=0.01
    #Add the following 2 lines.
    identity = float(sys.argv[3])
    pvalue = float(sys.argv[4])
    #-----20110930-------------------------------
    conDict = {}
    for line in open(sys.argv[1]):
        #to deal with long.table, some of them have no query sequence
        if line.strip() == '# BLAST processed 0 queries':
            #print >>sys.stderr, sys.argv[1], 'no long sequence'
            sys.exit(1)
        if line.startswith('AT'):
            line = line.strip()
            lineL = line.split('\t')
            query = lineL[0]
            hit = lineL[1]
            #patch for select identity and pvalue. 20110930
            #identity = 0.7 pvalue=0.01
            #Add the following 5 lines.
            if float(lineL[2]) < identity or \
                float(lineL[-2]) > pvalue:
                continue
            #Before add at 20110930
            if query not in conDict:
                conDict[query] = {}
                #print conDict
            if hit not in conDict[query]:
                conDict[query][hit] = set()
                #print conDict
            #print lineL[2:]
            conDict[query][hit].add('\t'.join(lineL[2:]))    
    #--------------------------------------
    #to deal with no hit, usually this will not happen.
    #if len(conDict) == 0:
    #    print >>sys.stderr, 'No hit', query
    #    return 0
    #rdict(conDict)
    subDict = {}
    for line in open(sys.argv[2]):
        if line[0] == '>':
            locus = line.rstrip()[1:]
            subDict[locus] = ''
        else:
            subDict[locus] += line.strip()
    numHomolog = len(subDict)
    #------------------------------------
    conservationRate(conDict, numHomolog-1)
if __name__ == '__main__':
    main()

