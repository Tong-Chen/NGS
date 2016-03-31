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
'''
Patch:
1.20110930--->add identity as a judge parameter in readMsaPsi and
readMsaPsiBatch. 
2.20110930--->add converge judgement. Before I assume all psiblast
are converged, then retun when find the word 'converged'. But it 
finally proved this is wrong, so if not converged, print to stderr,
and return when finish reading.
'''

import sys
import os

def readMsaPsi(file, identity=70, evalue=10):
    '''
    This function read the file produced by psiblast use MSA as the
    input. It returns a dict, query as a key, value is a set, '#'
    seperated subject, identity, start and end site of subject
    sequence.
    '''
    psiD = {}
    '''
    psiD = {locus, {sublocus:
        set("identity#start#end", "identity#start#end")}}
    '''
    filename = os.path.split(file)[1]
    if 0:
        print filename
    locus = '.'.join(filename.split('.', 3)[:-1])
    psiD[locus] = {}
    converged = 0
    for line in open(file):
        if line[0] != '#' and len(line) > 1:
            if line[:-1] == "Search has CONVERGED!":
                converged = 1
                return psiD
            else:
                lineL = line.split('\t')
                #print lineL
                if float(lineL[10]) <= evalue and \
                        float(lineL[2]) >= identity:
                    #----templ modofication----
                    hitlocus = lineL[1]
                    if hitlocus.find('|') != -1:
                        sub = hitlocus.split('|')[1]
                    else:
                        sub = hitlocus
                    #---templp modification--20111011
                    if sub not in psiD[locus]:
                        psiD[locus][sub] = set()
                    psiD[locus][sub].add('#'.join((\
                        lineL[2],lineL[8],lineL[9])))
            #-----------END one line-------------------
        #---------------END one line-------------------
    #-----------------END a file----------------
    if not converged:
        print >>sys.stderr, 'Uncoverged', file
        return psiD
#---------------readMsaPsi------------------------------

def readMsaPsiBatch(path, identity=70, evalue=10):
    '''
    This is a batch reader of readMsaPsi. It accept a path, and return a
    Dict.
    '''
    psiD = {}
    for file in os.listdir(path):
        file = path + file
        tmpD = readMsaPsi(file, identity, evalue)
        lenpsiDO = len(psiD)
        lentmpD = len(tmpD)
        psiD.update(tmpD)
        assert (lenpsiDO + lentmpD == len(psiD))
    #--------end----------------
    return psiD
#-------------------------------------------------

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
                pos = ((line.split('_')[0]).split('.'))[-1]
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

#def getSeqFromDB(locus, db):
#    '''
#    This function extracts sequence from formeted DB use 'blastdbcmd'
#    in BLAST+ package. 
#    give one locus, please makesure single is not 0. It will 
#    rueturn a fast formated sequences.  
#    '''
#    pass
#
#def getBatchSeqFromDB(file, db):
#    '''
#    This function extracts sequence from formeted DB use 'blastdbcmd'
#    in BLAST+ package. It needs the file contains the locus of
#    sequences wanted (one locus one line) as default. If you want to
#    give one locus, please makesure single is not 0. It will 
#    rueturn a fast formated sequences.  
#    '''
#    aDict = {}
#    if single:
#        cmd = "blastdbcmd -entry " + file + " -db " + db
#        seqList = os.popen(cmd).readlines()
#        lens = len(seqList)
#        if lens > 0:
#            key = seqList[0].strip()
#            aDict 
#            for i in range(1,lens):
#                aDict[key] +=
#    for locus in open(file):


if __name__ == '__main__':
    pass
