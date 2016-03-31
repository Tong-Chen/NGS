#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
import sys
'''
This is used to merge the repetitions together.


Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
==============================================================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=============================================================================================

#the format of key is locus, the value is a list. The element in value(list) is
#a dict which key is position, value is repetition. 

def readFromFile(filename, aDict):
    #-------------save information-----------------------------------
    #aDict = {key:[{position:repetition,position:repetition},{}]}
    for line in open(filename):
        if line.find('>') == 0:
            outerKey = line[1:-1]
            aDict[outerKey] = []

        if line.find('#') != -1:
            innerDict = {}
            for i in line.split('#'):
                if i.find(':') != -1: #'' is the last element
                    keyValuePair = i.split(':')
                    innerDict[int(keyValuePair[1])] = keyValuePair[0]
            aDict[outerKey].append(innerDict)
def merge(aDict):
    '''
    aDict = {outerKey:[{pos:repetition},{pos:repetition}]}
    '''
    for key, valueList in aDict.items():
        needmerge = 1
        valueListLen= len(valueList)
        if valueListLen > 1:
            mergeSamePos(valueList, needmerge, valueListLen)
    #------End of merge------------------------------------
def mergeSamePos(valueList, needmerge, valueListLen):
    '''
    valueList = [{pos:sequence},{pos:sequence}]
    '''
    needmerge = 0
    for no1 in range(valueListLen-1):
        dict1 = valueList[no1]
        for pos1 in dict1.keys():
            for no2 in range(no1+1, valueListLen):
                dict2 = valueList[no2]
                if dict2.has_key(pos1):
                    needmerge = 1
                    #print 'Merge'
                    #print dict1
                    #print dict2
                    joinDictSamePos(dict1,dict2)
                    break #this breaks the innest for
            #--------for no2 in range(no1, valueListLen)----
            if needmerge:
                break
        #-------for pos1 in dict1.keys()------
        if needmerge:
            break
    #----for no1 in range(valueListLen-1)-----    
    #remove the dict which value is '#'
    for dict0 in valueList:
        if '#' in dict0.values():
            try:
                assert ''.join(dict0.values()) == \
                    '#'*len(dict0.values())
            except AssertionError:
                print valueList
                return
            valueList.remove(dict0)
    if needmerge:
        valueListLen = len(valueList)
        if valueListLen > 1:
            mergeSamePos(valueList, needmerge, valueListLen)
#-------------mergeSamePos----------------------------



def joinDictSamePos(dict1, dict2):
    '''
    It joins two dicts which have repetitions at the same postition,
    and only save the longer repetitions.
    All the element will be saved in dict1, the dict.value in dict2
    will be '#' 
    '''
    for key, value in dict2.items():
        if dict1.has_key(key):
            if len(value) > len(dict1[key]):
                dict1[key] = value
        else:
            dict1[key] = value
        dict2[key] = '#'
#----joinDictSamePos(dict1, dict2)-----------

def output(filehandle, aDict):
    #--------------------output------------------------------------------
    laststd = sys.stdout
    sys.stdout = filehandle
    keylist = aDict.keys() #a list of protein locus
    keylist.sort()
    for i in keylist:
        sys.stdout.write('>%s\n' % i)
        #print len(aDict[i]) #how many epetitions
        for innerDict in aDict[i]:
            innerKey = innerDict.keys()
            #innerKeyNum = [int(keystr) for keystr in innerKey]
            innerKey.sort()
            for key in innerKey:
                #key = str(key)
                sys.stdout.write('%s:%s#' % (innerDict[key], key))
            sys.stdout.write('\n')
    sys.stdout = laststd

def main():
    length = len(sys.argv)
    if length != 3 and length != 2 :
        print 'Using %s inputfile [outputfile]' % sys.argv[0]
        print '---%s Tair.result.rearrange Tair.result.merge\n' % sys.argv[0]
        sys.exit(0)

    aDict = {}
    readFromFile(sys.argv[1], aDict)
    merge(aDict)
    if length == 3:
        fh = open(sys.argv[2], 'w')
        output(fh, aDict)
        fh.close()
    elif length == 2:
        output(sys.stdout, aDict)

if __name__ == '__main__':
    main()
