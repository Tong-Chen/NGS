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
    print >>sys.stderr, "Print the result to file."
    print >>sys.stderr, "Get the group with Arabidopsis thaliana expression data support."
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s groupname' % sys.argv[0]
        sys.exit(0)
    #----------------------------------------------------------------
    expressData = '/home/CT/server/project/project1/ori.material/\
expression/AtGE.DEV.gcRMA.AT.L.S'

    locusList = []
    for i in open(expressData):
        locusList.append(i[:-1])
    #-----------------------------------------------------------------
    lineDict = {}
    dictKey = 0
    supportedKey = set()
    for line in open(sys.argv[1]):
        if line[0] == '>':
            dictKey += 1
            lineDict[dictKey] = ''
        elif line.find('| Dicot.Arabidopsis.thaliana') != -1:
            num, locus = tuple(line.split()[0:2])
            locus = locus.split('.')[0]
            if locus in locusList:
                supportedKey.add(dictKey)
                line = line.replace(num, num+'*')
            #--------------------------------------------    
        lineDict[dictKey] += line

    lensupprt = len(supportedKey)

    print >>sys.stderr, "** %d supported group.**" % lensupprt
    
    if 1:
        output = sys.argv[1] + '.gcRMA'
        fh = open(output, 'w')

        for key in supportedKey:
            print >>fh, lineDict[key]

if __name__ == '__main__':
    main()

