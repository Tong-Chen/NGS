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
import os
from blast import readMsaPsiBatch

def getLocus(psiD, subjL):
    '''
    '''
    locusL = set()
    for valueD in psiD.values():
        [locusL.add(locus) for locus in valueD.keys()]
    #----------------------------------------------
    fh = open(subjL, 'w')
    print >>fh, '\n'.join(locusL),
    fh.close()

def output(psiD):
    for key, valueD in psiD.items():
        print '>%s' % key
        for vkey, vvalueS in valueD.items():
            print vkey, vvalueS

def readMsa(path):
    '''
    '''
    msaDict = {}
    for file in os.listdir(path):
        key = '.'.join(file.split('.', 3)[:-1])
        file = path + file
        tmpL = []
        tmpline = ''
        pos = ''
        for line in open(file):
            if line[0] == '>':
                if tmpline:
                    tmpline += ':' + str(pos)
                    tmpL.append(tmpline)
                tmpline = ''
                ##for in sativa there is '_' in locus. So changed the
                ##following line ---20111011
                #pos = ((line.split('_')[0]).split('.'))[-1]
                pos = ((line.rsplit('_', 1)[0]).split('.'))[-1]
                ##--20111011----
            else:
                #add 20110930, beefore forgot substitute -
                #patched 20111006, to deal with multiple lines.
                tmpline += line.strip().replace('-','')
            #-------------------------------
        #--------------------------------
        #--add the following 3 lines at 20111006. 
        #To deal with multiple lines. 
        if tmpline:
            tmpline += ':' + str(pos)
            tmpL.append(tmpline)
        #----------------------------------------
        tmpL.sort(key=lambda s: int(s.split(':')[1]))
        if key not in msaDict:
            msaDict[key] = '#'.join(tmpL)
        else:
            print >>sys.stderr, "Duplicated locus", key
        #----------End one file---------------
    #-------------All file in path------------------
    return msaDict
#-------------------------------------------------------------                
def readSubjS(subjS):
    '''
    '''
    subjSDict = {}
    for line in open(subjS):
        #print 'here'
        if line[0] == '>':
            ##--temp modification -20111011--
            #assert(line[:3] == ">gi")
            if line.find('|') != -1:
                gi = line.split('|')[1]
            else:
                gi = line.strip()[1:]
            ##--temp modification -20111011--
            if gi not in subjSDict:
                subjSDict[gi] = ''
            else:
                print >>sys.stderr, 'Duplicate gi', gi, 'in', subjS
        else:
            subjSDict[gi] += line.strip()
        #------End one line--------------------
    #----------End one file-------------------
    return subjSDict
#-------------------------------------------------

def deleteOverlap(alist, lenlist):
    '''
    alist = ['acgt:11', 'agcta':22]
    '''
    for i in range(0, lenlist-1):
        iseq, istart = alist[i].split(':')
        iend = int(istart) + len(iseq) - 1
        for j in range(i+1, lenlist):
            jseq, jstart = alist[j].split(':')
            jstart = int(jstart)
            if iend >= jstart:
                alist[i] = iseq + jseq[iend-jstart+1:] + ':' + istart
                alist[j] = alist[i]
            else:
                break
        #----------inner for---------------
    #--------------outer for--------------
    alist = set(alist)
    alist = list(alist)
    alist.sort(key=lambda s: int(s.split(':')[1]))
    ####patch---delete the same position one--
    tmpdict = {}
    for item in alist:
        seq, pos = item.split(':')
        if pos in tmpdict:
            if len(seq) > len(tmpdict[pos]):
                tmpdict[pos] = seq
        else:
            tmpdict[pos] = seq
    #---Before add 20111009----
    keyL = tmpdict.keys()
    keyL.sort(key=lambda s: int(s))
    alist = [':'.join((tmpdict[item], item)) for item in keyL]
    return alist
#---------deleteOverlap------------------------------
def deleteSamePosition(alist):
    tmpdict = {}
    for item in alist:
        seq, pos = item.split(':')
        if pos in tmpdict:
            if len(seq) > len(tmpdict[pos]):
                tmpdict[pos] = seq
        else:
            tmpdict[pos] = seq
    #---Before add 20111009----
    keyL = tmpdict.keys()
    keyL.sort(key=lambda s: int(s))
    alist = [':'.join((tmpdict[item], item)) for item in keyL]
    return alist


def getSeq(subjL, subjS,\
        db=r"/data7T/mercu-b-backup/pub-data/NCBI-NR/nr"):
    cmd = "blastdbcmd -entry_batch "+subjL+" -db " +db+' > '+ subjS
    os.system(cmd)

def final(psiD, subjSDict, msaDict, psiresult, delete):
    '''
    '''
    fh = open(psiresult, 'w')
    psiDKL = psiD.keys()
    psiDKL.sort()
    for pkey in psiDKL:
        #patched 20111006, add judgement of significant hit.
        pvalueD = psiD[pkey]
        #print >>fh, '=%s %d' % (pkey, len(pvalueD)) 
        #print >>fh, '>%s' % pkey
        rep142 = msaDict[pkey]
        #print >>fh, rep142
        minLen = 0.7 * len((rep142.split(':')[0]))
        pvalueDKL = pvalueD.keys()
        pvalueDKL.sort()
        pvalueDTmp = {}
        for pvkey in pvalueDKL:
            ####------
            seq = subjSDict[pvkey]
            pvalueDVL = list(pvalueD[pvkey])
            pvalueDVL.sort(key=lambda s: int(s.split('#')[1]))
            pvtmpL = []
            for pvitem in pvalueDVL:
                iden,start,end = pvitem.split('#')
                start = int(start)
                end = int(end)
                if (end-start+1) >= minLen:
                    rep = seq[start-1:end] + ':' + str(start)
                    pvtmpL.append(rep)
            #--------End extract one subject-----
            ###The following 4 lines added 20111008. To deal
            ###with the overlapped hits. 
            if delete:
                lenpvtmpl = len(pvtmpL)
                if lenpvtmpl > 1:
                    pvtmpL = deleteOverlap(pvtmpL, lenpvtmpl)
            ###Lines before and function add 20111008.
            ###add 20111010 delete same position but different length
            else:
                pvtmpL = deleteSamePosition(pvtmpL)
            ###line before changed 20111010
            ##lines after changed 20111009----to delete short hits--
            if len(pvtmpL):
                pvalueDTmp[pvkey] = '#'.join(pvtmpL)
            #print >>fh, '>%s' % pvkey
            #print >>fh, '#'.join(pvtmpL)
            #-------END one subbject---
        #-----------END all subject---
        lenpvalueDTmp = len(pvalueDTmp)
        print >>fh, '=%s %d' % (pkey, lenpvalueDTmp) 
        print >>fh, '>%s' % pkey
        print >>fh, rep142
        pvalueDTmpKeyL = pvalueDTmp.keys()
        if lenpvalueDTmp:
            ##for not all locus like gi are all numbers, so no int
            #pvalueDTmpKeyL.sort(key=lambda s: int(s))
            pvalueDTmpKeyL.sort()
            ##for not all locus like gi are all numbers, so no int
            ##above changed 20111011
            for item in pvalueDTmpKeyL:
                print >>fh, '>%s\n%s' % (item, pvalueDTmp[item])
        #-----------END one key----
    #---------------END all
    fh.close()
#---------------------------------------------

def main():
    print >>sys.stderr, "Print the result to file"
    if len(sys.argv) < 8:
        print >>sys.stderr, 'Using python %s psipath/ psiresult \
subjectlocus subjectSeq msaPath/ identity[70] evalue[10] \
deleteOverlap[0-nodelete, 1-delete(default)]{The last parameter is \
optional}' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------------------------------
    psipath = sys.argv[1]
    psiresult = sys.argv[2]
    subjL = sys.argv[3]
    subjS = sys.argv[4]
    msaPath = sys.argv[5]
    identity = float(sys.argv[6])
    evalue = float(sys.argv[7])
    if len(sys.argv) == 9:
        delete = int(sys.argv[8])
    else:
        delete = 1
    psiD = readMsaPsiBatch(psipath, identity, evalue)
    #output(psiD)
    #------------------------------------------------
    #pathced 20111006, add judgement of existing files
    nolocus = 1
    if os.path.isfile(subjL) and os.path.isfile(subjS):
        nolocus = 0
        nolocus = int(raw_input("Do you want to recreate hitted locus\
            file and seq file? Default no[0], 1 means yes. \n>>>"))
    if nolocus:
        getLocus(psiD, subjL)
        getSeq(subjL, subjS)
    #-----------------------------------------
    #print subjS
    subjSDict = readSubjS(subjS)
    #print subjSDict
    #sys.exit(1)
    msaDict = readMsa(msaPath)
    #print subjSDict
    final(psiD, subjSDict, msaDict, psiresult, delete)
if __name__ == '__main__':
    main()

