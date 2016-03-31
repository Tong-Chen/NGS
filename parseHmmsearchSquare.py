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
def readHmm(hmmout):
    fh = open(hmmout)
    resutDict = {}
    while 1:
        line = fh.readline()
        if len(line) == 0: break
        while line[0:3] != ">> ":
            line = fh.readline()
            if len(line) == 0: break
        if len(line) == 0: break
        locus = line.split()[1]
        line = fh.readline() #   #      score  bias****************
        if len(line) == 0: break
        line = fh.readline() # ---     --------*
        if len(line) == 0: break
        line = fh.readline() #the first line of result
        if len(line) == 0: break
        while len(line) > 1:
            #-------------save hit-----------
            if line.find('[]') != -1:
                if locus not in resutDict:
                    resutDict[locus] = 1
                else:
                    resutDict[locus] += 1
            #-------------save hit-----------
            line = fh.readline()
            if len(line) == 0: break
        #------one hit------------------------------
    #----------total hit--------------------------
    fh.close()
    return resutDict
#----------------End readHmm---------------------------
def main():
    print >>sys.stderr, "This is used to extract the hitted proteins \
of a motif using hmmsearch.Here we only take the full length match \
as positive result. Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #----------------------------------------------------
    resutDict = readHmm(sys.argv[1])
    keyL = resutDict.keys()
    keyL.sort()
    print '\n'.join(keyL)
    print >>sys.stderr, "Total hits %s proteins %s hits" % \
            (len(keyL), sum(resutDict.values()))

if __name__ == '__main__':
    main()

