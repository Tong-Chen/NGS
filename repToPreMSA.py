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
#import os
from ctIO import readRep
def main():
    print >>sys.stderr, "Transfer repetitions to multiple sequence \
alignment files. One geoup in one file. Waiting for alignment."
    if len(sys.argv) != 2:
        print >>sys.stderr,'Using python %s repfile' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------
    repDict = {}
    readRep(sys.argv[1], repDict)
    for locus, valueL in repDict.items():
        midlen = 30
        group = 0
        for groupD in valueL:
            group += 1
            groupDKeyL = groupD.keys()
            groupDKeyL.sort()
            maxlen = 0
            for seq in groupD.values():
                lenseq = len(seq)
                if lenseq > maxlen:
                    maxlen = lenseq
            if maxlen <= midlen:
                file = locus+'.'+str(group)+'.short'
            else:
                file = locus+'.'+str(group)+'.long'
            fh = open(file, 'w')
            for pos in groupDKeyL:
                seq = groupD[pos]
            #-------------------------------------------
                posn = ':'.join((str(pos[0]), str(pos[1])))
                print >>fh, '>%s.%s.%s\n%s' % \
                    (locus, str(group), posn, seq)
            #--------END one group------------------------------
            fh.close()
#            msashort = 'clustalw ' + short + ' -OUTPUT=FASTA'
#            msalong = 'clustalw ' + long + ' -OUTPUT=FASTA'
#            os.system(msashort)
#            os.system(msalong)
        #------------END one locus
    #-------------END----all-----------------
if __name__ == '__main__':
    main()

