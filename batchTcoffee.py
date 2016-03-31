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
import os
from ctIO import readRep
def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s repfile' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------------
    repdict = {}
    readRep(sys.argv[1], repdict)
    for id, valueL in repdict.items():
        i = 0
        for seqDict in valueL:
            filename = id + str(i) + '.fasta'
            fh = open(filename, 'w')
            seqDictK = seqDict.keys()
            seqDictK.sort()
            output = ''
            for key in seqDictK:
                output += ''.join(('>pos', str(key[0]), '-', str(key[1])\
                    , '\n', seqDict[key], '\n'))
            #--------------------------------------------
            print >>fh, output,
            fh.close()
            cmd = 't_coffee ' + filename
            i += 1
            os.system(cmd)
        #---------------End one id-------------- 
    #--------End all id-------------------------

if __name__ == '__main__':
    main()

