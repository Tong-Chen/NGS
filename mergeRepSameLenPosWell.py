#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
This is used to merge the repetitions found to form same length
sequences.

Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
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
            ####print line,

        if line.find('#') != -1:
            innerDict = {}
            keyValuePair = []
            ####saveI = []
            for i in line[:-1].rstrip('#').split('#'):
                if i.find(':') == -1:
                    #Some accident, some repeat lost the position
                    #information
                    print >>sys.stderr, 'wrong'
                    return
                    #i = getIndex(outerKey, seqDict, i, keyValuePair)
                ####saveI.append(i)
                keyValuePair = i.split(':') #keyValuePair = [pep, pos]
                innerDict[int(keyValuePair[1])] = keyValuePair[0]
            aDict[outerKey].append(innerDict)
            del innerDict
            ####print '#'.join(saveI)

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
                if int(num) >= small:
                    return i+':'+str(num)
    print indexes, small
    return 'wrong'
#--------------end-------------------------------------
def merge(aDict):
    '''
    aDict = {outerkey:[{pos:sequence}, {pos:sequence}]}
    '''
    for key, valueList in aDict.items():
        needmerge = 1
        #print 'Before Processing >%s\n@%s' % (key, valueList)
        valueListLen = len(valueList)
        if valueListLen > 1:
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
        lendict1 = len(dict1)
        for no2 in range(no1+1, valueListLen):
            dict2 = valueList[no2]
            if lendict1 == len(dict2):
                pos1arr = dict1.keys()
                pos1arr.sort()
                pos2arr = dict2.keys()
                pos2arr.sort()
                for posi in range(lendict1):
                    needmerge = \
                        posCompare(pos1arr[posi],pos2arr[posi],
                            len(dict1[pos1arr[posi]]),\
                            len(dict2[pos2arr[posi]]))
                    if not needmerge: #if there is one element which
                                    #not need merge, break 
                        break
                #---end for about comparing whether it can merge
                if needmerge:
                    #print 'Before merge'
                    #print dict1
                    #print dict2
                    joinDict(dict1, dict2, pos1arr, pos2arr, lendict1)
                    #print 'After merge'
                    #print dict1
                    #print dict2
                    break #get out of the second iteration
        #---------end of second iteration----for ----        
        if needmerge:
            break  #get out of the first iteration
    #----end of first iteration ---for----------
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
            dealMerge(valueList, needmerge, valueListLen)
#----------------------end-------------------------------------

def joinDict(dict1, dict2, pos1arr, pos2arr, lendict1):
    '''
    joinDict(dict1, dict2, pos1arr, pos2arr, lendict1)

    this joins two dicts together if these two dicts have get 
    through the test.

    It only joins two related sequences together.
    All the element will be saved in dict1, the values of dict2 will
    be '#'.
    '''
    for posi in range(lendict1):
        joinPep(pos1arr[posi],pos2arr[posi],dict1, dict2,\
                dict1[pos1arr[posi]], dict2[pos2arr[posi]],\
                len(dict1[pos1arr[posi]]), len(dict2[pos2arr[posi]]))

def posCompare(pos1, pos2, len1, len2):
    '''
    posCompare(pos1, pos2, len1, len2)-->int(0 or 1)

    this used to decide whether this two seq can merge,
    '''
    if pos1 <= pos2 < pos1+len1:
        return 1
    elif pos2 <= pos1 < pos2+len2:
        return 1
    return 0


def joinPep(pos1, pos2, dict1, dict2, pep1, pep2, len1, len2):
    '''
    joinPep(pos1, pos2, dict1, dict2, pep1, pep2, len1, len2)
        
    Join two sequences.
    '''
    if pos1 <= pos2 < pos1+len1:
        if pos2+len2 > pos1+len1:
            dict1[pos1] = \
                    pep1[0:pos2-pos1]+pep2
        dict2[pos2] = '#' 
    elif pos2 <= pos1 < pos2+len2:
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
            innerKey.sort()
            for key in innerKey:
                sys.stdout.write('%s:%s#' % (innerDict[key], key))
            sys.stdout.write('\n')
    sys.stdout = laststd

def main():
    print >>sys.stderr, 'This is used to merge repetitions, it gets new\
        repetitions in same length.'
    length = len(sys.argv)
    if length != 2 and length != 3 :
        print 'Using %s repfile [outputfile]' % sys.argv[0]
        print '>>>%s TAIR9.RESULT.50.4 Tair.result.merge\n' \
                % sys.argv[0]
        sys.exit(0)

    #the format of key is locus, the value is a list. The element in value(list) is
    #a dict which key is position, value is the length of repetition. 
    aDict = {}
    #seqDict = {}
    
    #readSeq(sys.argv[2], seqDict)
    readFromFile(sys.argv[1], aDict)
    merge(aDict)
    if length == 3:
        fh = open(sys.argv[2], 'w')
        output(fh, aDict)
        fh.close()
    elif length == 2:
        output(sys.stdout, aDict)
        pass

if __name__ == '__main__':
    main()
