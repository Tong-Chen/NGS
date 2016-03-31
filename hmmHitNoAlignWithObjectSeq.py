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

def readObjectSeq(file):
    '''
    #file format
    >locus
    ATCGATCGCAGCAT
    aDict = {locus:seq}
    '''
    aDict = {}
    for line in open(file):
        if line[0] == '>':
            locus = line[1:-1]
            assert (locus not in aDict)
        else:
            aDict[locus] = line.strip()
    #---------------------------------------
    return aDict
#---------------------------------------
def readHmm(hmmfile, include=0.1):
    '''
    one hmmhit one file

    incluse significant hits and hits with inclusion threshold less
    than 0.1(or other include value given) 
    and its group has already contained. No length information
    needed for the using of di-orientation varification.
    '''
    aDict = {}
    #next = 1000000
    #new = 0
    fh = open(hmmfile)
    line = fh.readline()
    while(not line.startswith("Query")):
        line = fh.readline()
    #--------curline strats with Query---
    print 'Query', line,
    nouse, searchLoc, searchLLen = line.strip().split()
    assert (searchLoc not in aDict)
    hitDict = {}
    aDict[searchLoc] = hitDict
    #hitDict = {{hitgroup:[hitloc, hitloc2, ]}, {}}
    i = 4
    while i:
        line = fh.readline()
        i -= 1
    #----------------------------
    line = fh.readline()
    print 'significant', line
    while line.find("inclusion") == -1:
        hitLoc = line.strip().rsplit(' ', 1)[1]
        hitLocGroup = hitLoc.split('_')[0]
        if hitLocGroup in hitDict:
            hitDict[hitLocGroup].append(hitLoc)
        else:
            hitDict[hitLocGroup] = [hitLoc]
        #-----------------------------------
        line = fh.readline()
        while len(line) == 1: #break when meet blank line and read next
            #line, which contains the inclusion  
            line = fh.readline()
    #--END significant--------------------------------
    print 'end significant', line
    #---begin inclusion---
    line = fh.readline()
    while len(line) == 1:
        line = fh.readline()
    print 'begin inclusion', line
    while not line.startswith("Domain annotation for each"):
        #print line
        #----link into--if one in this group hitted, take
        #others in if in threshold with evalue less 
        #than 0.1 
        hitLoc = line.strip().rsplit(' ', 1)[1]
        hitLocGroup = hitLoc.split('_')[0]
        if hitLocGroup in hitDict:
            evalue = float(line[0:12]) #emperical number
            if evalue < include: #inclusion evalue
                hitDict[hitLocGroup].append(hitLoc)
        #-----------------------------
        line = fh.readline()
        while len(line) == 1: 
            #break when meet blank line and read next
            #line, which contains the inclusion  
            line = fh.readline()
    #-------------End inclusion--------------------
    print 'End inclusion', line
    return aDict
#---------------------------------------------------


def readHmmSimple(hmmfile):
    '''
    one hmmhit one file

    incluse significant hits and hits with inclusion threshold less
    than 0.1(or other include value given) 
    and its group has already contained. No length information
    needed for the using of di-orientation varification.
    '''
    aDict = {}
    #next = 1000000
    #new = 0
    fh = open(hmmfile)
    line = fh.readline()
    while(not line.startswith("Query")):
        line = fh.readline()
    #--------curline strats with Query---
#    print 'Query', line,
    nouse, searchLoc, searchLLen = line.strip().split()
    assert (searchLoc not in aDict)
    #hitDict = {}
    aDict[searchLoc] = set()
    #hitDict = {{hitgroup:[hitloc, hitloc2, ]}, {}}
    i = 4
    while i:
        line = fh.readline()
        i -= 1
    #----------------------------
    line = fh.readline()
    while len(line) == 1:
        line = fh.readline()
#    print 'significant', line
    while line.find("inclusion") == -1 and \
        line.find('No hits') == -1:
        #print line, hmmfile
        hitLoc = line.strip().rsplit(' ', 1)[1]
        hitLocGroup = hitLoc.split('_')[0]
        aDict[searchLoc].add(hitLocGroup)
        #if hitLocGroup in hitDict:
        #    hitDict[hitLocGroup].append(hitLoc)
        #else:
        #    hitDict[hitLocGroup] = [hitLoc]
        #-----------------------------------
        line = fh.readline()
        if len(line) == 1: #break when meet blank line and read next
            #line, which contains the inclusion  
            #line = fh.readline()
            break
    #--END significant--------------------------------
#    print 'end significant', line
#    #---begin inclusion---
#    line = fh.readline()
#    while len(line) == 1:
#        line = fh.readline()
#    print 'begin inclusion', line
#    while not line.startswith("Domain annotation for each"):
#        #print line
#        #----link into--if one in this group hitted, take
#        #others in if in threshold with evalue less 
#        #than 0.1 
#        hitLoc = line.strip().rsplit(' ', 1)[1]
#        hitLocGroup = hitLoc.split('_')[0]
#        if hitLocGroup in hitDict:
#            evalue = float(line[0:12]) #emperical number
#            if evalue < include: #inclusion evalue
#                hitDict[hitLocGroup].append(hitLoc)
#        #-----------------------------
#        line = fh.readline()
#        while len(line) == 1: 
#            #break when meet blank line and read next
#            #line, which contains the inclusion  
#            line = fh.readline()
#    #-------------End inclusion--------------------
#    print 'End inclusion', line
    fh.close()
    return aDict
#---------------------------------------------------

def secAnalysisSimple(aDict):
    '''
    aDict={searchLoc:(hitgroup,) }
    '''
    bDict = {}
    for searchLoc, hitSet in aDict.items():
        for hitgroup in hitSet:
            if hitgroup != searchLoc and hitgroup in aDict:
                if searchLoc in aDict[hitgroup]:
                    if searchLoc not in bDict:
                        bDict[searchLoc] = set()
                    bDict[searchLoc].add(hitgroup)
    #-------------------------------------------------
    return bDict
#--------------End secAnalysisSimple------------------
def thirdSimple(bDict):
    '''
    bDict = {locus, set(hit1, hit2,)}
    '''
    cDict = {}
    for locus, valueS in bDict.items():
        already = 0
        #-------judge if this locus already in key or value
        if len(cDict):
            for locusS, valueSS in cDict.items():
                if locus in valueSS:
                    already = 1
                    break
        #---------add new key-----------------
        if not already:
            cDict[locus] = set()
            for item in valueS:
                linkItem(item, bDict, locus, cDict)
        #--------------End one locus---------------------
    #---------END for---------------------
    return cDict
#-------------END function---------------

def linkItem(item, bDict, locus, cDict):
    ''''''
    if item in cDict[locus] or item == locus:
        return
    cDict[locus].add(item)
    itemValueS = bDict[item]
    for itemValueSItem in itemValueS:
        if itemValueSItem in cDict[locus] or \
            itemValueSItem == locus:
            continue
        else:
            cDict[locus].add(itemValueSItem) 
            for itemValueSItemSItem in bDict[itemValueSItem]:
                linkItem(itemValueSItemSItem, bDict, locus, cDict)
        #------------------------------------            
#-----------------END function--------

def outputSimple(aDict):
    ''''''
    for item, valueS in aDict.items():
        itemg = item.split('.',1)[0]
        label=''
        aset = set()
        aset.add(itemg)
        for valueSI in valueS:
            valueSIPre = valueSI.split('.', 1)[0]
            aset.add(valueSIPre)
        #---------------------------------------
        len247 = len(aset)
        if len247 > 1:
            label = '*'
        #------------------------------------
        print ">%s%s %d" % (item, label, len247) 
        print '#'.join(valueS)
#---------outputSimple----------    

#def output(aDict):
#    '''
#    hitDict = {searchLoc : {
#        hitgroup:[hitloc, hitloc2, ], 
#        hitgroup:[hitloc, hitloc2, ]
#          }
#    '''
#    newDict = {}
#    for  searchLoc, hitDict in aDict.items():
#        for hitDictK, hitDictVL in hitDict.items():
#            if hitDictK != searchLoc and hitDictK in aDict:
#                hitDictKValueD = aDict[hitDictK]
#                if searchLoc in hitDictKValueD.keys():
#                    for i in 
##----------------------------------

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s hmmhitpath/' % sys.argv[0]
        sys.exit(0)
    #----------------------------
    import os
    hmmDIct = {}
    path = sys.argv[1]
    for file in os.listdir(path):
        if file.endswith('Hit'):
            file = path + file
            aDict = readHmmSimple(file)
            hmmDIct.update(aDict)
    #--------------------------------------
    #print hmmDIct
    hmmDIctSec = secAnalysisSimple(hmmDIct)
    hmmDIctThird = thirdSimple(hmmDIctSec)
    outputSimple(hmmDIctThird)
if __name__ == '__main__':
    main()

