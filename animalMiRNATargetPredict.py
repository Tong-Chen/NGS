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
from time import localtime, strftime 
import string
timeformat = "%Y-%m-%d %H:%M:%S"

def readtarget(target):
    targetDict = {}
    for line in open(target):
        if line[0] == '>':
            locus = line[1:-1]
            assert locus not in targetDict
            targetDict[locus] = ''
        else:
            targetDict[locus] += line.strip().upper().replace('T','U')
    #-------------------------------------------
    return targetDict
#-------------------------------------------------

def alignOneMirnaMultipleTarget(targetDict, seed, hitNum):
    target_pre = []
    for key, seq in targetDict.items():
        if seq.count(seed) >= hitNum:
            target_pre.append(key)    
    #----------------------------------------
    return target_pre
#----------------------------------------

def readMir(mirna):
    #------------------------------
    #>cel-miR-1-5p MIMAT0020301 Caenorhabditis elegans miR-1-5p
    #CAUACUUCCUUACAUGCCCAUA
    #------------------------------
    mirnaDict = {}
    for line in open(mirna):
        if line[0] == '>':
            locus = line[1:].strip().split()[0]
            mirnaDict[locus] = ''
        else:
            mirnaDict[locus] += line.strip().upper()
    #-------------------------------
    return mirnaDict
#--------------------------------------------
def main():
    '''
    1.character case, all transfer to uppercase
    2.all U and T will transfer to A in mirna seed reverse
    complementary.
    3.all utr sequence T->U
    '''
    if len(sys.argv) != 5:
        print >>sys.stderr, "This is used to predict miRNA target by \
seed regions. If given miRNA has at least <hitNum> matchs of seed \
region in given fasta seq, this will be taken as a miRNA-gene pair."
        print >>sys.stderr, 'Using python %s matureMirna[fasta,\
sequence in one line] \
targetSequence[fasta, 3"UTR, CDS, 5"UTR] hitNum outputPrefix' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    transTable = string.maketrans('ACGUT', 'UGCAA')
    mirna=sys.argv[1]
    target=sys.argv[2]
    hitNum=int(sys.argv[3])
    targetDict = readtarget(target)
    mirnaDict = readMir(mirna)
    #---mirna index----------------------
    newmirnaDict = {} #for target index
    file = sys.argv[4]+ '.' + str(hitNum) + '.miRNA.index.target' 
    fh = open(file, 'w')
    for key, value in mirnaDict.items():
        value = value[1:8]
        value = value.translate(transTable)
        value = value[::-1]
        newmirnaDict[key] = value
        target_pre = alignOneMirnaMultipleTarget\
            (targetDict, value,hitNum)
        for i in target_pre:
            print >>fh, '%s\t%s' % (key, i)
    fh.close()
    #-------targetIndex--------------------------
    file = sys.argv[4] + '.target.index.mirna' 
    fh = open(file, 'w')
    for key, value in targetDict.items():
        mirna_pre = []
        for mirna_key, mirna_value in newmirnaDict.items():
            count88 = value.count(mirna_value)
            if count88:
                mirna_pre.append(mirna_key+'\t'+str(count88))
            #-------------------------------------
        #-------------------------------------
        for item in mirna_pre:
            print >>fh, "%s\t%s" % (key, item)
    fh.close()
    
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    #startTime = localtime()
    #startTime = '-'.join([str(x) for x in startTime[:3]]) \
    #    + ' ' + ':'.join([str(x) for x in startTime[3:6]])
    main()
    endTime = strftime(timeformat, localtime())
    #endTime = localtime()
    #endTime = '-'.join([str(x) for x in endTime[:3]]) \
    #    + ' ' + ':'.join([str(x) for x in endTime[3:6]])
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


