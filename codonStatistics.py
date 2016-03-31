#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
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
    '''
    Detect the codon usage in different repeats.
    '''
    print >>sys.stderr, "Print the result to files"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    dataDict = {}
    file=sys.argv[1]
    for line in open(file):
        if line[0] == '>':
            pass
            #print line,
        else:
            replist = [seg.split(':')[0] for seg in \
                    line.rstrip('\n#').split('#')]
            lenL = len(replist[0])
            assert lenL % 3 == 0
            lenL = int(lenL / 3)
            count = 0  #record number of different codon at same
                       #position
            for i in range(lenL):
                tmplist = []
                tmplist = [item[i*3:(i+1)*3] for item in replist]
                tmplist = set(tmplist)
                if len(tmplist) > 1:
                    count += 1
                #------------------
            #---------------------
            #print count, lenL
            #assert isinstance(count, int)
            #assert isinstance(lenL, int)
            #print '%f' % (count / lenL)
            diff = count / lenL
            if diff not in dataDict:
                dataDict[diff] = 1
            else:
                dataDict[diff] += 1
        #----End one line-------------------------
    #--------End reading--------------------------
    filehist=file +'.Hist'
    fh = open(filehist, 'w')
    for key, value in dataDict.items():
        for time in range(value):
            print >>fh, key
    #---------------------------------
    fh.close()
    #---------------------------------
    fileaccum = file + '.Accum'
    fh = open(fileaccum, 'w')
    accumuDict = {}
    accumuDictKl = []
    anniDIct = {}
    step = 0.05
    sum = int(1/step + 1.2) #1.2 just for in case, 1 is ok 
    init = 0
    for i in range(sum):
        accumuDict[init] = 0
        anniDIct[init] = 0
        accumuDictKl.append(init)
        init += step
    #--------------------------
    for key, value in dataDict.items():
        for i in accumuDictKl:
            if key <= i:
                accumuDict[i] += value
    #----------------------------------
    #----animation data---------------
    subtraction = accumuDict[accumuDictKl[-1]] - \
        accumuDict[accumuDictKl[0]]
    idealstep = subtraction / (sum -1)
    anniDIct[accumuDictKl[0]] = accumuDict[accumuDictKl[0]]
    anniDIct[accumuDictKl[-1]] = accumuDict[accumuDictKl[-1]]
    for i in range(1, sum-1):
        anniDIct[accumuDictKl[i]] = anniDIct[accumuDictKl[i-1]]\
            + idealstep 
    #----animation data---------------
    for i in accumuDictKl:
        print >>fh, "%.2f %d %.2f" % (i, accumuDict[i], anniDIct[i])
    fh.close()
if __name__ == '__main__':
    main()

