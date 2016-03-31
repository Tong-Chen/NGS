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
from ctIO import readFasta, readRep
#from ctTEST import ct_rdict

def codonSet():
    str = "ATCG"
    codonList = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codonList.append(str[i]+str[j]+str[k])
    #-------------------------------------
    codonList.sort()
    return codonList
#----------------------------------------------------

def originalSta(repDict, cdsDict):
    codonRepDictTmp = {}
    '''codonRepDictTmp={'locus': {'triplestartpos':codon}}'''
    codonSeqDict = {}
    '''codonSeqDict={'locus': {'codon':number}}'''
    for locus, seq in cdsDict.items():
        codonSeqDict[locus] = {}
        lenseq = len(seq)
        #A special locus--------------
        if lenseq %3 != 0:
            lenseq -= lenseq%3
            print locus
            print seq
        assert lenseq % 3 == 0
        for i in range(0, lenseq, 3):
            tmpkey1 = seq[i:i+3]
            if tmpkey1 in codonSeqDict[locus]:
                codonSeqDict[locus][tmpkey1] += 1
            else:
                codonSeqDict[locus][tmpkey1] = 1
        #--End seq----------------------------------------
        #patch a bug, before repDict part has indented 4 more spaces
        #forget why I did that before. 2011-08-25
        if locus in repDict:
            codonRepDictTmp[locus] = {}
            for valueD in repDict[locus]:
                for valueDKey in valueD.keys():
                #the number in valueDKey is the real position
                    start = valueDKey[0]-1
                    end = valueDKey[1]-1
                    for i in range(start, end, 3):
                        if i not in codonRepDictTmp[locus]:
                            codonRepDictTmp[locus][i] =\
                                seq[i:i+3]
                    #-------------------------------------
        #----------End---------rep------------------------------
    #--------------End---------seq----------------------------
    #------------------------finish---------------------------
    #------------------------transfer----------------------
    codonRepDict = {} #'''codonRepDict={codon:num}'''
    for keytmp, dicttmp in codonRepDictTmp.items():
       codonRepDict[keytmp] = {}
       for dicttmpvalue in dicttmp.values():
           if dicttmpvalue not in codonRepDict[keytmp]:
               codonRepDict[keytmp][dicttmpvalue] = 1
           else:
               codonRepDict[keytmp][dicttmpvalue] += 1
    #-----------------finish----transfer----------------------
    return codonRepDict, codonSeqDict

#--------------------------------------------------------------
def sta(srcDict, objDict):
    for codon, num in srcDict.items():
        if codon not in objDict:
            objDict[codon] = num
        else:
            objDict[codon] += num
#------------------------------------------------
def singlRatioProRep(codonRepDict, codonSeqDict, file, codonSet):
    '''
    Divide codon num of repeat segment by codon num of other segement.
    Output a matrix, row is 64 codons in each protein, colum is each
    codon in all proteins. 
    '''
#    codonSet = set()
#    for valueD in codonSeqDict.values():
#        [codonSet.add(codon) for codon in valueD.keys()]
#    codonSet = list(codonSet)
#    codonSet.sort()
    file = file + '.Heatmap'
    fh = open(file, 'w')
    print >>fh, '\t'.join(codonSet)
    for locus, valueD in codonRepDict.items():
        tmpList = []
        sumup = sum(valueD.values())
        sumdown = sum(codonSeqDict[locus].values())
        sumdown -= sumup
        for codon in codonSet:
            up = valueD[codon] if codon in valueD else 0
            down = codonSeqDict[locus][codon] if codon in \
                codonSeqDict[locus] else 0  # here down can be any
                #number except 0, for down is denominator, up is 0
                #I changed, there is a situation that codon only
                #exists in 
                #ratio
            #------this codon does not exist-----
            if down == 0:
                if up == 0:
                    tmpList.append('NA')
                    continue
                else:
                    print >>sys.stderr, "Wrong"
                    sys.exit(1)
            #------this codon does not exist-----
            down -= up
            #---deal with diff length------------
            up /= sumup
            down /= sumdown
            #---if codon only exists in rep---------
            if down == 0 and up != 0:
                tmpList.append('Inf')
                continue
            #if down == 0:
            #    print up, down, valueD, codon, codonSeqDict[locus]
            tmpList.append(str(up/down))
        print >>fh, "%s\t%s" % (locus, '\t'.join(tmpList))
    fh.close()

def totalNumberProrep(codonRepDict, codonSeqDict, file):
    '''
    This computs the total number of 64 codons in repeat
    coding segements and other parts of proteins which have
    repeats.
    '''
    #patch to exclude stop codons 2011-08-25
    stopCodon = ['TAA', 'TAG', 'TGA']
    codonNumRep = {}
    codonNumSeq = {}
    
    for locus, valueD in codonRepDict.items():
        sta(valueD, codonNumRep)
        sta(codonSeqDict[locus], codonNumSeq)
    #--------------------------------------------
    for codon, num in codonNumSeq.items():
        diff = codonNumRep[codon] if codon in codonNumRep else 0
        codonNumSeq[codon] -= diff
    #-------------------------------------------
    #print '#codon rep seq'
    file = file +  '.CodonRepSeq.BiBar'
    fh = open(file, 'w')
    print >>fh, "# rep other"
    codonSet = codonNumSeq.keys()
    codonSet.sort()
    numrep = sum(codonNumRep.values())
    numseq = sum(codonNumSeq.values())
    for codon in codonSet:
        if codon in codonNumRep:
            rep = codonNumRep[codon]
        else:
            rep = 0
        seq = codonNumSeq[codon]
        #patch to exclude stop codons 2011-08-25
        if codon not in stopCodon:
            print >>fh, '%s %f %f' % (codon, rep/numrep, seq/numseq)
    fh.close()
    return codonNumSeq
#------------------------------------------------------------

def repOrNot(codonNumSeq, codonRepDict, codonSeqDict, file, codonList):
    '''
    This compares codon usage between repeat-containg proteins(exclude
    repeat part) and repeat-no proteins.
    
    A bug: total length different. So output ratio of one codon of 64
    codons in one group.
    '''
    #patch to exclude stop codons 2011-08-25
    stopCodon = ['TAA', 'TAG', 'TGA']
    repList = codonRepDict.keys()
    #this saves the proteins do not containg rep
    codonNumSeqNorep = {}

    for locus, valueD in codonSeqDict.items():
        if locus not in repList:
            sta(valueD, codonNumSeqNorep)
    #------------------------------------------------------
    file = file +  'Codon.RepNo.Bibar'
    fh = open(file, 'w')
    print >>fh, "# rep norep"
    sumrep = sum(codonNumSeq.values())
    sumno = sum(codonNumSeqNorep.values())
    for codon in codonList:
        rep = codonNumSeq[codon] if codon in codonNumSeq else 0
        no = codonNumSeqNorep[codon] if codon in codonNumSeqNorep else 0
        #patch to exclude stop codons 2011-08-25
        if codon not in stopCodon:
            print >>fh, "%s %f %f" % (codon, rep/sumrep, no/sumno)
    fh.close()
#------------------------------------------------------
def main():
    print >>sys.stderr, "Print the result to three files"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s seq rep' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------------
    codonList = codonSet()
    cdsDict = readFasta(sys.argv[1])
    #ct_rdict(cdsDict)
    repDict = {}
    readRep(sys.argv[2], repDict)
    #ct_rdict(repDict)
    codonRepDict, codonSeqDict = originalSta(repDict, cdsDict)
    #ct_rdict(codonRepDict)
    #print '*********************'
    #ct_rdict(codonSeqDict)
    #-------compare within protein with repeats----
    #--get codons within repeat and divide codons within other
    #seuquences, and bar graph the number of them, heatmap the 
    #ratio of each codons of one protein.
    codonNumSeq = \
        totalNumberProrep(codonRepDict, codonSeqDict, sys.argv[2])
    singlRatioProRep(codonRepDict, codonSeqDict, sys.argv[2], codonList)
    #--compare proteins have no repeat and proteins have repeats but
    #delete repeats  
    repOrNot(codonNumSeq, codonRepDict, codonSeqDict, sys.argv[2],
            codonList)

    
if __name__ == '__main__':
    main()

