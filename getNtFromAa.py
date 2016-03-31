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

#--------patch a bug---2011-08-25

import sys
import ctIO

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s cdsfile repfile' \
            % sys.argv[0]
        sys.exit(0)
    #---------------------------------------------------------
    repDict = {}
    ctIO.readRep(sys.argv[2], repDict)
    locusL = repDict.keys()
    locusL.sort()
    cdsDict = ctIO.readFasta(sys.argv[1], locusL)

    for locus in locusL:
        print '>%s' % locus
        seq = cdsDict[locus]
        tmpList = repDict[locus]
        for posDict in tmpList:
            posKeys = posDict.keys()
            posKeys.sort()
            repList = []
            for posTuple in posKeys:
                start = (posTuple[0] - 1) * 3
                end = posTuple[1] * 3
                if start >= end:
                    print >>sys.stderr, locus, posTuple
                    sys.exit(1)
                #--------patch a bug---2011-08-25
                #repList.append(seq[start:end]+':'+str(start+3))
                repList.append(seq[start:end]+':'+str(start+1))
            #--------------------------------------------------
            print '#'.join(repList)
        #------------------------------------------------------
    #----------------------------------------------------------

if __name__ == '__main__':
    main()

