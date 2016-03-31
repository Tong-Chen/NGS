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

def alignOneMirnaMultipleTarget(mirna, targetDict, mature, transTable):
    seed = mature[1:8]
    seed = seed.translate(transTable)
    seed = seed[::-1]
    lenm = len(mature)
    newmature = mature[::-1]
    for key, seq in targetDict.items():
        target_pre = []
        index = seq.find(seed)
        lenseq = len(seq)
        while index != -1:
            utr_s = index - lenm + 8 
            utr_e = index + 8 #not included
            #------------------
            newpos_s = index
            newpos_e = index + 7
#            if utr_s < 0:
#                newpos_s = 0
#            else:
#                newpos_s = utr_s
#            if utr_e > lenseq:
#                newpos_e = lenseq
#            else:
#                newpos_e = utr_e
            matchSeq = "Position " + str(newpos_s+1) + '-' + str(newpos_e)\
                + ' of ' + key + ' 3\' UTR'
            #------------------
            matchSeq += "\n\n5\'  ..."
            if utr_s < 0:
                utr_s = 0
                matchSeq += '-' * (-1 * utr_s)
            #-------------------------
            if utr_e > lenseq:
                seq += '-' * (lenseq - utr_e)
            matchSeq += seq[utr_s:utr_e]
            matchSeq += '\n       ' + ' ' * (lenm-8) + '|' * len(seed) 
            matchSeq += '\n3\'     ' + newmature
            target_pre.append(matchSeq)    
            oldindex = index
            index = seq[index+1:].find(seed)
            if index != -1:
                index += oldindex+1
        #------------END while--------------
        if len(target_pre) > 0:
            print "%s\n\n%s\n\n" % (mirna, '\n\n'.join(target_pre))
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
    print >>sys.stderr, "Print the result to file"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s matureMirna[fasta,\
sequence in one line] \
targetSequence[fasta, 3"UTR, CDS, 5"UTR]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    transTable = string.maketrans('ACGUT', 'UGCAA')
    mirna=sys.argv[1]
    target=sys.argv[2]
    targetDict = readtarget(target)
    mirnaDict = readMir(mirna)
    #---mirna index----------------------
    #newmirnaDict = {} #for target index
    for key, value in mirnaDict.items():
        #newmirnaDict[key] = value
        alignOneMirnaMultipleTarget\
            (key, targetDict, value, transTable)
    
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


