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

#from ctIO import readFasta

def readMeme(memeOutput):
    aDict = {}
    fh = open(memeOutput)
    label1 = 'BL   MOTIF '
    label2 = '//'
    while 1:
        line = fh.readline()
        if len(line) == 0: break
        #-----start-----------------------
        while not line.startswith(label1):
            line = fh.readline()
            if len(line) == 0: break
        if len(line) == 0: break
        key = line.split()[2]
        if key not in aDict:
            aDict[key] = {}
        else:
            print "Duplicated key %s" % key
            sys.exit(1)
        line = fh.readline()
        if len(line) == 0: break
        #----end and extract-------------------
        while not line.startswith(label2):
            locus, null, null, seq, null = line.split() 
            if locus not in aDict[key]:
                aDict[key][locus] = seq
            else:
                print "Duplicated locus %s" % locus
                sys.exit(1)
            line = fh.readline()
        if len(line) == 0: break
        #-----------------------------------------------------
    #------------------------------------------------------
    fh.close()
    return aDict
#----------------------------------------------------

def outputMemeDict(memeDict, outputdir):
    for key, valueD in memeDict.items():
        fh = open(outputdir+'motif_'+key, 'w')
        locusKeys = valueD.keys()
        locusKeys.sort()
        for locus in locusKeys:
            print >>fh, '>%s\n%s' % (locus, valueD[locus])
        fh.close()
#------------------------------------------

def main():
    print >>sys.stderr, "Get each motif from meme.txt, and save it \
in the file named after motif_num in FASTA format. Here the file is \
assumed to be the clustralw fasta format output."
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s meme.txt outputdir/' % sys.argv[0]
        sys.exit(0)
    #--------------------------------------------------------
    #memeInput = sys.argv[1]
    memeOutput = sys.argv[1]
    outputdir = sys.argv[2]
    #repDict = readFasta(memeInput)
    memeDict = readMeme(memeOutput)
    outputMemeDict(memeDict, outputdir)
if __name__ == '__main__':
    main()

