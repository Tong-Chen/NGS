#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division #, with_statement
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
    print >>sys.stderr, "Summarize average coverage of each feature, \
usually input file is the output from coverageBed and a gene feature."
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 4:
        print >>sys.stderr, 'Using python %s filename width[height] \
prefix' % sys.argv[0]
        sys.exit(0)
    #------------------------------------
    if sys.argv[2] == 'width':
        covIndex = -3
    elif sys.argv[2] == 'height':
        covIndex = -4
    #---------------------------------------
    aDict ={}
    '''
    aDict = {'gene':{genename:[cov, width]}, 'I':{genename:[cov, width]}}
    '''
    for line in open(sys.argv[1]):
        lineL = line.split()
        key = lineL[3]
        if key.find('.') != -1:
            fkey, skey = key.split('.', 1)
        else:
            fkey = key
            skey = 'gene'
        cov = int(lineL[covIndex])
        width = int(lineL[-2])
        if skey not in aDict:
            aDict[skey] = {}
        if fkey not in aDict[skey]:
            aDict[skey][fkey] = [cov, width]
        #---------------------------------------
        else:
            if skey == 'I' or skey == 'E':
                aDict[skey][fkey][0] += cov
                aDict[skey][fkey][1] += width
            else:
                print >>sys.stderr, "Wrong ",  fkey,  skey
                sys.exit()
        #----------------------------------------
    #-----------------------------------------
    #----------output--------------------
    prefix = sys.argv[3]
    skeyL = aDict.keys()
    skeyL.sort()
    for item in skeyL:
        fh = open(prefix+'.'+item, 'w')
        fkeyL = aDict[item].keys()
        fkeyL.sort()
        for geneName in fkeyL:
            print >>fh,  "%s\t%f" % (geneName, \
                aDict[item][geneName][0]/aDict[item][geneName][1])
        fh.close()
    #-----------------------------------------------------------
if __name__ == '__main__':
    main()

