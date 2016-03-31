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
from ctIO import readRep, outputRep
from ctAssist import shannonIndex as si

def saveDict(aDict, key, value):
    if key not in aDict:
        aDict[key] = [value]
    else:
        aDict[key].append(value)
#----------------------------------------
#def output(aDict, file):
#    locusL = aDict.keys()
#    locusL.sort()
#    fh = open(file, 'w')
#    for key in locusL:
#        print >>fh, '>%s' % key
#        itemDL = aDict[key]
#        for itemD in itemDL:
#            groupL = []
#            itemDKeyL = itemD.keys()
#            itemDKeyL.sort()
#            for posset in itemDKeyL:
#                groupL.append(':'.join((itemD[posset],str(posset[0]))))
#            #-------------------------------------------
#            print >>fh, '#'.join(groupL)
#        #-----------------------------------------------------
#    #---------------------------------------------------------
#    fh.close()
##----------------------------------------
def main():
    print >>sys.stderr, "Using the average shannonIndex value \
of a group sequences to represent the last entropy."
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s filename\
 threshold(2)[threethe more the high complexity]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    #this three dict have the same structure
    repDict = {}
    lcsDict = {} #save low complexity sequences
    regularDict = {} #save regular sequences
    readRep(sys.argv[1], repDict)
    lcs = int(sys.argv[2])
    for locus, valueL in repDict.items():
        for itemD in valueL:
            entropy = 0
            i_valueS = set(itemD.values())
            #i_keys = itemD.keys()
            for item in i_valueS:
                entropy += si(item)
            entropy = entropy / len(i_valueS)
            if entropy <= lcs:
                saveDict(lcsDict, locus, itemD)
            else:
                saveDict(regularDict, locus, itemD)
        #--------End one dict---------------
    #-------------end all-----------------
    prefile = sys.argv[1].split('/')[-1]
    outputRep(lcsDict, prefile+'.LCSs')
    outputRep(regularDict, prefile+'.HCSs')
if __name__ == '__main__':
    main()

