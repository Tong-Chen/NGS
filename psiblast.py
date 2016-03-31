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

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------------------
    hitDict = {}
    tmpline = ''
    for line in open(sys.argv[1]):
        if line.startswith('Query= '):
            locus = line.strip().split()[1]
            #print '***%s\n\n' % locus
            hitDict[locus] = {}
        elif line[0] == '>':
            tmpline = line.strip()
        elif line.startswith('Length='):
            if tmpline:
                #print tmpline
                tmplineL = tmpline.split('|')[1:]
                lentmplineL = len(tmplineL)
                print tmplineL
                for i in range(0, lentmplineL, 2):
                    hitlocus = tmplineL[i]
                    taxonmyTmp = tmplineL[i+1]
                    if taxonmyTmp.find('[') != -1 and \
                        taxonmyTmp.find(']') != -1:
                        hittaxonmy = \
                            (taxonmyTmp.split('[')[1]).split(']')[0]
                    else:
                        hittaxonmy = 'unknown'
                    if hittaxonmy in hitDict[locus]:
                        hitDict[locus][hittaxonmy].append(hitlocus)
                    else:
                        hitDict[locus][hittaxonmy]=[]
            #-----------------------------------------
            tmpline = ''
        elif tmpline:
            tmpline += line.strip()
    #--------------------------------------------------------------
    #print '*' * 50
    #print hitDict
    #print '*' * 50
    for locus, adict in hitDict.items():
        keyL = adict.keys()
        keyL.sort()
        print locus, keyL, len(keyL)
if __name__ == '__main__':
    main()

