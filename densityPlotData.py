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

def main():
    print >>sys.stderr, "Statistics the repetition density along the protein sequence"
    if len(sys.argv) != 3:
        print 'Using python %s FOR.C merge' % sys.argv[0]
        sys.exit(0)
    
    locusLenDict = {}

    flag = 0
    locus = ''
    for line in open(sys.argv[1]):
        if flag == 0:
            locus = line[:-1]
        elif flag == 1:
            locusLenDict[locus] = int(line[:-1])
        elif flag == 2:
            flag = -1
        flag += 1
    #---------end for-------------------------
    
    locusDistriDict = {}
    distriList = []
    for line in open(sys.argv[2]):
        if line[0] == '>':
            distriList = []
            locus = line[1:-1]
            locusDistriDict[locus] = []
        else:
            distriList = [(int(seqpos.split(':')[1]),\
                int(seqpos.split(':')[1])+ len(seqpos.split(':')[0]))\
                for seqpos in line[:-1].rstrip('#').split('#')]
            locusDistriDict[locus].append(distriList)

    #-------------end for------------------------
    #locusDensity = {}
    for locus, distriListList in locusDistriDict.items():
        length = locusLenDict[locus]
        density = [0] * length
        for distriList in distriListList:
            for posRange in distriList:
                for pos in range(posRange[0], posRange[1]):
                    density[pos-1] = density[pos-1] + 1
                    #print '%d\t%d\n' % (pos, density[pos-1])
                #--------------for---------------------------
            #------------------for---------------------------
        #----------------------for---------------------------
        #locusDensity[locus] = density
        filename = locus
        fh = open(filename, 'w')
        #j = 1
        for i in density:
            print >>fh, i
            #j += 1
        fh.close()
    #--------------------------for---------------------------
    
if __name__ == '__main__':
    main()

