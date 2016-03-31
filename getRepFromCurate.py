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
import ctIO

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s seqforc curatedRep' % sys.argv[0]
        sys.exit(0)
    seqDict = {}
    ctIO.readseq(sys.argv[1], seqDict)
    repDict = {}
    for line in open(sys.argv[2]):
        if line[0] == '>':
            locus = line[1:]
            locus = locus.strip()
            seq = seqDict[locus]
            shift = 0
            repDict[locus] = {}
        else:
            starPos = line.find('*')
            if starPos != -1:
                rep = line[:starPos]
                rep = rep.strip()
                rep = rep.replace(' ','')
                repstart = seq.find(rep, shift)
                if repstart == -1:
                    print >>sys.stderr, 'No such repetition %s'\
                        % rep
                    sys.exit(0)
                shift = repstart + len(rep)
                key = '*' * line.count('*')
                if key not in repDict[locus]:
                    repDict[locus][key] = []
                repDict[locus][key].append(rep+':'+str(repstart+1))
            #---------id rep-------------------
        #----------Not locus line-------------
    #-------------End reading-----------------
    for locus, sonDict in repDict.items():
        if len(sonDict) > 0:
            print '>%s' % locus
            for value in sonDict.values():
                print '#'.join(value)
        else:
            print >>sys.stderr, locus

            
if __name__ == '__main__':
    main()

