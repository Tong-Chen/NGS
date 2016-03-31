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

def test(aDict, bDict, core):
    '''
    dict = {'A':('B1','B2',...)}
    '''
#    for aItem, valueS in aDict.items():
#        for bItem in valueS:
#            if aItem not in bDict[bItem]:
#                print >>sys.stderr, aItem, bItem
    if aDict[core] != bDict[core]:
        print >>sys.stderr, aDict[core]
        print >>sys.stderr, bDict[core]
        print >>sys.stderr, 'Asymmetric'
    for itemB in aDict[core]:
        if itemB == core:
            continue
        itemAS = bDict[itemB]
        if len(itemAS) != 1 or core not in itemAS:
            print >>sys.stderr, itemAS, itemB
#---------------------------------------------
def readBiogrid(file, col, core='', head =1):
    '''
    core is the molecular you used to get the interaction file. If no
    core is given, do not do any test.
    col = [1,2],[3,4],[5,6],[7,8](Starts begin 0, 
    not the real column number)
    '''
    aDict = {}
    bDict = {}
    for line in open(file):
        if head:
            head -= 1
            continue
        interactorA, interactorB = \
            line.split('\t')[col[0]:col[1]+1]
        #--------------biline-----------------
        #If these two rows are added, no use the sentence in fucntion
        #test.
        if interactorA == interactorB:
            continue
        #--------------biline-----------------
        if interactorA not in aDict:
            aDict[interactorA] = set()
        aDict[interactorA].add(interactorB)
        if interactorB not in bDict:
            bDict[interactorB] = set()
        bDict[interactorB].add(interactorA)
        #-----------------------------------------
    #--------------------------------------------
    #if core:
    #    test(aDict, bDict, core)
    #    test(bDict, aDict, core)
    #----------------------------------------
    newset = aDict[core]
    newset.update(bDict[core])
    return newset
#-------------------------------------------------
def outputResult(aDict, core):
    ''''''
    print '\n'.join(aDict[core])

#-------------------------------------
def main():
    print >>sys.stderr, "This is used to extract the interaction \
proteins of your given [core] protein from a biogrid2 format file \
downloaded as searches result form the server."
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 4:
        print >>sys.stderr, 'Using python %s filename core test[0 and\
 other strings mean "false" for no test]' % sys.argv[0]
        sys.exit(0)
    #--------------------------------------------
    core = sys.argv[2]
    test = sys.argv[3]
    col = [5,6]
    if test:
        print '\n'.join(readBiogrid(sys.argv[1], col, core))
    else:
        print '\n'.join(readBiogrid(sys.argv[1], col))
    #-----------------------------------------------

#--------------------------------------------------------------
if __name__ == '__main__':
    main()

