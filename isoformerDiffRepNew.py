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
from ctIO import readIsoformId, readAnno

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) > 4:
        print >>sys.stderr, 'Using python %s repId locusHasIsoform \
            [annofile]' % sys.argv[0]
        sys.exit(0)
    #----------------------------------------------------
    allIsoformId = readIsoformId(sys.argv[2])
    repIsoformId = readIsoformId(sys.argv[1])
    hasAnno = 0
    if len(sys.argv) == 4:
        hasAnno = 1
        annodict = {}
        readAnno(sys.argv[3], annodict)

    repKL = repIsoformId.keys()
    repKL.sort()
    print "WithRep\tNoRep\tAnno\tGene Ontology\tInterpro"
    for repK in repKL:
        if repK not in allIsoformId:
            continue
        repVL = repIsoformId[repK]
        lenRepVL = len(repVL)
        allVL = allIsoformId[repK]
        lenAllVL = len(allVL)
        #----------diff-------------------------
        if lenRepVL < lenAllVL:
            repVL.sort()
            allVL.sort()
            #------output-----------------------
            for allV in allVL:
                locus = repK+allV
                print >>sys.stderr, locus
                if hasAnno:
                    if locus in annodict:
                        anno = annodict[locus]
                        anno = anno.replace(r'\\', '\t')
                    if allV in repVL:
                        print "%s\t*\t%s" % (locus, anno)
                    else:
                        print "*\t%s\t%s" % (locus, anno)
                else:
                    if allV in repVL:
                        print "%s\t*" % locus
                    else:
                        print "*\t%s" % locus

            #------output-----------------------
        #----------diff-------------------------
    #--------------trace-----------------------

if __name__ == '__main__':
    main()

