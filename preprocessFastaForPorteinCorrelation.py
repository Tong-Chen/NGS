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
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"

def readFasta(file):
    key = ''
    aDict = {}
    for line in open(file):
        if line[0] == '>':
            if key:
                aDict[key] = ''.join(aDict[key])
            key=line.strip()[1:]
            if key not in aDict:
                aDict[key] = []
            else:
                print >>sys.stderr, "Duplicated key %s" % key
        else:
            aDict[key].append(line.strip())
        #------------------------------------------------------------
    if key:
        aDict[key] = ''.join(aDict[key])
    #---end reading---------------------------------------------------------
    #-----------compare----------------------------------------------------
    bDict = {}
    for key, value in aDict.items():
        tmpkey = key.split(':',1)[0]
        if tmpkey not in bDict:
            bDict[tmpkey] = value
        else:
            if len(value) > len(bDict[tmpkey]):
                bDict[tmpkey] = value
    #-----------------------------------------------
    return bDict
#-----------------------END-------------------------------------

def delete_same_func(aDict, lenAL, aveA, bDict, lenBL, aveB, comKey, lenKeyCom, delete_same):
    delL = []
    for i in range(lenKeyCom - 1):
        if i in delL:
            continue
        keyi = comKey[i]
        val_i = aDict[keyi]
        for j in range(i+1, lenKeyCom):
            if j in delL:
                continue
            keyj = comKey[j]
            val_j = bDict[keyj]
            #----same on set A---------------------
            if val_i == val_j:
                if delete_same == '2':
                    #---------same also in set B-------
                    if bDict[keyi] == bDict[keyj]:
                        delL.append(j)
                elif delete_same == '1':
                    #---choose the one which length close to average
                    #length----
                    if abs(lenBL[i]-aveB)>abs(lenBL[j]-aveB):
                        delL.append(i)
                        break
                    else:
                        delL.append(j)
            #---------------------------------------
        #-------inner for-------------------------
    #-----------outter for------------------------
    #----------delete bDict if delete_same == '1'------
    if delete_same == '1':
        for i in range(lenKeyCom - 1):
            if i in delL:
                continue
            keyi = comKey[i]
            val_i = bDict[keyi]
            for j in range(i, lenKeyCom):
                if j in delL:
                    continue
                keyj = comKey[j]
                val_j = aDict[keyj]
                #----same on set A---------------------
                if val_i == val_j:
                        #---choose the one which length close to average
                        #length----
                        if abs(lenAL[i]-aveA)>abs(lenAL[j]-aveA):
                            delL.append(i)
                            break
                #---------------------------------------
            #-------inner for-------------------------
        #-----------outter for------------------------
    #---------------END delete ------------------------       
    #--------update aDict, bDict------------------------
    for i in delL:
        aDict.pop(comKey[i])
        bDict.pop(comKey[i])
    #-----------------------------------------------------

#----------------------END delete_same_func----------------------

def main():
    lensysargv = len(sys.argv)
    if lensysargv != 3 and lensysargv != 4:
        print >>sys.stderr, 'Using python %s fasta1 \
fasta2[>species:otherthings, usually orthoDB output file] \
delete_same[default 2. Accept 0 or 2. 0 means not delete same \
sequences. 1 means delete sequences same in one dataset. 2 means \
delete sequences same in both datasets.] \
' % sys.argv[0]
        sys.exit(0)
    #---unused--------------------------------
    delete_same = '2'
    if lensysargv > 3:
        delete_same = sys.argv[3]
    #---unused--------------------------------
    aDict = readFasta(sys.argv[1])
    bDict = readFasta(sys.argv[2])
    #----------same species---------------
    keyA = set(aDict.keys())
    keyB = set(bDict.keys())
    comKey = list(keyA.intersection(keyB))
    comKey.sort()
    for key in keyA:
        if key not in comKey:
            aDict.pop(key)
    for key in keyB:
        if key not in comKey:
            bDict.pop(key)
    #---------get average length---------
    lenKeyCom = len(comKey)
    if lenKeyCom < 2:
        print >>sys.stderr, "Less sequences"
        sys.exit(1)
    lenAL = [len(value) for value in aDict.values()]
    aveA  = sum(lenAL) / lenKeyCom
    lenBL = [len(value) for value in bDict.values()]
    aveB  = sum(lenBL) / lenKeyCom
    #---------delete same----------------
    if delete_same:
        delete_same_func(aDict, lenAL, aveA, bDict, lenBL, aveB, comKey, lenKeyCom, delete_same)
        keyA = set(aDict.keys())
        keyB = set(bDict.keys())
        comKey = list(keyA.intersection(keyB))
        comKey.sort()
    #----------output----------------------
    output1 = sys.argv[1]+'_'+sys.argv[2]+'.Rev'
    fh = open(output1, 'w')
    for key in comKey:
        print >>fh, ">%s\n%s" % (key, aDict[key])
    fh.close()
    output2 = sys.argv[2]+'_'+sys.argv[1]+'.Rev'
    fh = open(output2, 'w')
    for key in comKey:
        print >>fh, ">%s\n%s" % (key, bDict[key])
    fh.close()
#----------------------------------------------------


if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


