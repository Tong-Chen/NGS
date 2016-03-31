#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''


Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
#from __future__ import division, with_statement
import sys
from ctIO import readRep

def transfer(repDict):
    '''
    original:(Fuzzy or exact)
    >AT1G01490.1
    GPAKEPEKEKKEE:69#PPKKEGEAPKEEG:90#KEGEAPKKEEEKK:104#KEGGDKKEGEKK:116#
    PKKEGEAPKEEG:91#KKEGEAPKKEEE:103#KKEGGDKKEGEK:115#
    transfered:
    >locus.grounum.multipleposition(: separate)
    >AT1G01490.1.1.69:90  #if 69 and 90 have the same sequence.
    GPAKEPEKEKKEE
    >AT1G01490.1.1.90
    PPKKEGEAPKEEG
    ...
    >AT1G01490.1.2.91
    PKKEGEAPKEEG
    '''
    repDictKeyL = repDict.keys()
    repDictKeyL.sort()
    for locus in repDictKeyL:
        groupL = repDict[locus]
        group = 0
        for groupD in groupL:
            group += 1
            tmpdict = {}
            groupDKeyL = groupD.keys()
            groupDKeyL.sort()
            for pos in groupDKeyL:
                seq = groupD[pos]
                if seq not in tmpdict:
                    tmpdict[seq] = [str(pos[0])]
                else:
                    tmpdict[seq].append(str(pos[0]))
            #---------------------------------------------
            tmpdictKeyL = tmpdict.keys()
            tmpdictKeyL.sort()
            for seq in tmpdictKeyL:
                pos = ':'.join(tmpdict[seq])
                print '>%s.%s.%s\n%s' % \
                    (locus, str(group), pos, seq)
            #---END one seq------------------------------------------
        #-------END one group---------------------
    #-----------END all------------------------
#-------------------------------------------------
def main():
    print >>sys.stderr, 'used to transfer merged repetition file \
to fasta, the output is STDOUT.'
    if len(sys.argv) != 2:
        print 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    repDict = {}
    readRep(sys.argv[1], repDict)
    transfer(repDict)
    
if __name__ == '__main__':
    main()

