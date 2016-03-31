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

def readSeq(filename, seqDict):
    '''
    AT1G51370.2
    346
    MVGGKKKTKICDKVSHEEDRISQLPEPLISEILFHLSTKDSVRTS
    '''
    isSeq = 0
    for line in open(filename):
        if isSeq == 0 and line[:-1].isalnum():
            isSeq = 1
        elif not isSeq:
            locus = line[:-1]
        elif isSeq:
            seqDict[locus] = line[:-1]
            isSeq = 0

def readFromFile(filename, aDict):
    #-------------save information-----------------------------------
    '''
    readFromFile(filename, aDict)
    
    aDict = {outerkey:[{pos:sequence}, {pos:sequence}]}
    '''
    for line in open(filename):
        if line.find('>') == 0:
            outerKey = line[1:-1]
            aDict[outerKey] = []

        if line.find('#') != -1:
            innerDict = {}
            keyValuePair = []
            for i in line[:-1].rstrip('#').split('#'):
                if i.find(':') == -1:
                    #Some accident, some repeat lost the position
                    #information
                    print 'wrong'
                    sys.exit(1)
                    #i = getIndex(outerKey, seqDict, i, keyValuePair)
                keyValuePair = i.split(':') #keyValuePair = [pep, pos]
                innerDict[int(keyValuePair[1])] = keyValuePair[0]
            aDict[outerKey].append(innerDict)
            del innerDict

def getIndex(outerkey, seqDict, i, keyValuePair):
    '''
    
    
    '''
    import kmp
    indexes = kmp.kmp_matcher(seqDict[outerkey], i)
    if indexes != -1:
        if len(indexes) == 1:
            return i+':'+str(indexes[0])
        else:
            small = int(keyValuePair[1])+len(keyValuePair[0])
            for num in indexes:
                if num > small:
                    return i+':'+str(num)
#--------------end-------------------------------------
def merge(aDict, choose):
    '''
    aDict = {outerkey:[{pos:sequence}, {pos:sequence}]}
    '''
    for key,  valueList in aDict.items():
        needmerge = 1
        #print 'Before Processing >%s\n@%s' % (key, valueList)
        valueListLen = len(valueList)
        if valueListLen > 1:
            if choose == '1':
                import mergeRepSameLen
                mergeRepSameLen.mergeSameLen(valueList, needmerge,\
                        valueListLen)
            elif choose == '2':
                import mergeSamePos
                mergeSamePos.mergeSamePos(valueList, needmerge, \
                        valueListLen)
            elif choose == '3':
                dealMerge(valueList, needmerge, valueListLen)
        #print 'After Processing >%s\n@%s' % (key, valueList)

def dealMerge(valueList, needmerge, valueListLen):
    '''
    dealMerge(valueList)

    valueList = [{pos:sequence}, {pos:sequence}]

    this merged the neighbored repeat together.

    Attention: this will lose some variation information.

    The merge conditions are listed as follow:neighbored sequence,
    part-overlap sequence, total-contain sequence, the same sequence.

    such as 
    (1)neighbored sequence  -> n
    ------------
                -----
    -----------------
    indexSmall < indexLarge == indexSmall+len(pepSmall)
    (2)part-overlap sequence  -> p
        ---------
    ---------
    or
    -----
    -----------
    indexSmall < indexLarge < indexSmall+len(pepSmall) <= indexLarge+len(pepLarge)
    or
    indexSmall == indexLarge and len(pepSmall) != len(pepLarge)
    (3)total-contain sequence  -> t
         ------
    -----------------
    indexSmall < indexLarge < indexLarge+len(pepSmall) < indexSmall+len(pepLarge)
    (4)the same end sequence  -> s
    -----------------
    -----------------
    or
    ------
    --------------
    or
           -------
    --------------
    indexSmall == indexLarge 

    The scheme is that:

    compare the index of the sequence, if the index satisfy one of the last 4
    conditions, then merge. When the index in dict1 is over, do the
    same operation for dict1.
    '''
    #if not needmerge:
    #    return
    needmerge = 0
    for no1 in range(valueListLen -1):
        dict1 = valueList[no1]
        for pos1, pep1 in dict1.items():
            assert isinstance(pos1, int)
            assert isinstance(pep1, str)
            #this dict has been deleted
            if pep1 == '#':
                break
            len1 = len(pep1)
            for no2 in range(no1+1, valueListLen):
                dict2 = valueList[no2]
                for pos2, pep2 in dict2.items():
                    #this dict has been deleted
                    if pep2 == '#':
                        break
                    len2 = len(pep2)
                    needmerge = posCompare(pos1,pos2,dict1,dict2,\
                            pep1, pep2, len1, len2)
                    if needmerge:
                        joinDict(dict1, dict2)
                        break
                #---------a for---------------------------
                if needmerge:
                    break
            #--------a for----------------------------
            if needmerge:
                break
        #---------a for-----------
        if needmerge:
            break
    #-----a for---------------
    for dict0 in valueList:
        if '#' in dict0.values():
            assert ''.join(dict0.values()) == \
                '#'*len(dict0.values())
            valueList.remove(dict0)
    if needmerge:
        valueListLen = len(valueList)
        if valueListLen > 1:
            dealMerge(valueList, needmerge, valueListLen)
#----------------------end-------------------------------------

def joinDict(dict1, dict2):
    '''
    joinDict(dict1, dict2)

    this joins two dicts together if these two dicts have common
    element.
    The join can complete in three steps:
        first, merge those that can.(this can be ignore)
        second,move those that can not.
        third, merge those in fict1.
    All the element will be saved in dict1, the values of dict2 will
    be '#'.
    '''
    for pos2, pep2 in dict2.items():
        if pep2 != '#':
            if pos2 in dict1:
                if len(pep2) > len(dict1[pos2]):
                    dict1[pos2] = pep2
            else:
                dict1[pos2] = pep2
            dict2[pos2] = '#'
    selfMerge(dict1)

def selfMerge(dict1):
    '''
    selfMerge(dict1)

    This does the self merging.
    '''
    posList = dict1.keys()
    posList.sort()
    lenPosList = len(posList)
    needmerge = 0
    for no1 in range(lenPosList-1):
        pos1 = posList[no1]
        pep1 = dict1[pos1]
        len1 = len(pep1)
        if pep1 == '#':
            continue
        for no2 in range(no1+1, lenPosList):
            pos2 = posList[no2]
            pep2 = dict1[pos2]
            len2 = len(pep2)
            if pep2 == '#':
                continue
            needmerge = posCompare(pos1, pos2, dict1, dict1, pep1,\
                    pep2, len1, len2)
            if needmerge:
                break
        #-----------inner for------------------
        if needmerge:
            break
    #---outer for----------------------
    if needmerge:
        selfMerge(dict1)
    else:
        for pos, pep in dict1.items():
            if pep == '#':
                del dict1[pos]
        return

def posCompare(pos1, pos2, dict1, dict2, pep1, pep2, len1, len2):
    '''
    posCompare(pos1, pos2, dict1, dict2, pep1, pep2, len1, len2)
        -->int(0 or 1)

    if we do not want neighbored sequence without overlap merged
    together, minus 1 of last element.
    '''
    needmerge = 0
    if pos1 <= pos2 < pos1+len1-1:
        needmerge = 1
        if pos2+len2 > pos1+len1:
            dict1[pos1] = \
                    pep1[0:pos2-pos1]+pep2
        dict2[pos2] = '#' 
    elif pos2 <= pos1 < pos2+len2-1:
        needmerge = 1
        del dict1[pos1]
        #pos1 = pos2
        if pos1+len1 > pos2+len2:
            newpep = pep2[0:pos1-pos2]+pep1
        elif pos1+len1 <= pos2+len2:
            newpep = pep2
        if pos2 in dict1:
            if len(dict1[pos2]) < len(newpep):
                dict1[pos2] = newpep
        else:
            dict1[pos2] = newpep
        dict2[pos2] = '#'
    return needmerge
#----------------end------------------


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
            #print >>sys.stderr, innerDict
            for key in innerKey:
                sys.stdout.write('%s:%s#' % (innerDict[key], key))
            sys.stdout.write('\n')
    sys.stdout = laststd

def main():
    print 'This is used to merge repetitions, but it gets new\
        repetitions in different length.'
    length = len(sys.argv)
    if length != 3 and length != 2 :
        print 'Using %s inputfile [outputfile]' % sys.argv[0]
        print '>>>%s Tair.result.rearrange Tair.result.merge\n' % sys.argv[0]
        sys.exit(0)
    print '''Please choose a number to do what types of merge.
        >>>1:mergeRepSameLen, 2:mergeSamePos, 4:mergeTotal,  
        >>>if you want to use two program you need add the number
        >>>together to get a sum as input'''
    choose = raw_input('Please input one number from 1,2,3,4,5,6\n>>>')
    #the format of key is locus, the value is a list. The element in value(list) is
    #a dict which key is position, value is the length of repetition. 
    aDict = {}
    #seqDict = {}
    
    #readSeq('TAIR8_pep_20080412_FOR_C', seqDict)
    readFromFile(sys.argv[1], aDict)
    merge(aDict, choose)
    if length == 3:
        fh = open(sys.argv[2], 'w')
        output(fh, aDict)
        fh.close()
    elif length == 2:
        output(sys.stdout, aDict)

if __name__ == '__main__':
    main()
