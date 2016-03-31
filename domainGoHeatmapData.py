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
import re

TESTdict=0
if TESTdict:
    from ctTEST import ct_rdict

'''test str
str1='Molecular Function: small conjugating protein ligase activity
(GO:0019787), Biological Process: post-translational protein
modification (GO:0043687), Biological Process: regulation of protein
metabolic process (GO:0051246)'

'''
pattern = re.compile(r' ?([^,]+?):(.+?)\((GO:.+?)\)')
def testpattern():
    '''
>>> str1
'Molecular Function: small conjugating protein ligase activity (GO:0019787), Biological Process: post-translational protein modification (GO:0043687), Biological Process: regulation of protein metabolic process (GO:0051246)'
>>> pat.findall(str1)
[('Molecular Function', ' small conjugating protein ligase activity ', 'GO:0019787'), ('Biological Process', ' post-translational protein modification ', 'GO:0043687'), ('Biological Process', ' regulation of protein metabolic process ', 'GO:0051246')]
>>> str2
'Molecular Function: two-component response regulator activity (GO:0000156), Biological Process: two-component signal transduction system (phosphorelay) (GO:0000160), Biological Process: regulation of transcription, DNA-dependent (GO:0006355)'
>>> pat.findall(str2)
[('Molecular Function', ' two-component response regulator activity ', 'GO:0000156'), ('Biological Process', ' two-component signal transduction system (phosphorelay) ', 'GO:0000160'), ('Biological Process', ' regulation of transcription, DNA-dependent ', 'GO:0006355')]
    '''
    str1 = 'Molecular Function: small conjugating protein ligase activity (GO:0019787), Biological Process: post-translational protein modification (GO:0043687), Biological Process: regulation of protein metabolic process (GO:0051246)'
    print pattern.findall(str1)
    str2 = 'Molecular Function: two-component response regulator activity (GO:0000156), Biological Process: two-component signal transduction system (phosphorelay) (GO:0000160), Biological Process: regulation of transcription, DNA-dependent (GO:0006355)'
    print pattern.findall(str2)
#-----------------End test-------------------------------
def addgo(go, bpDict, ccDict, mfDict, locus):
    for item in pattern.findall(go):
        if item[0] == 'Molecular Function':
            mfdsp = item[1]
            mfdsp = mfdsp.strip()
            mfsym = item[2]
            add2dict2(mfDict, mfsym, mfdsp, locus)
        elif item[0] == 'Biological Process':
            bpdsp = item[1]
            bpdsp = bpdsp.strip()
            bpsym = item[2]
            add2dict2(bpDict, bpsym, bpdsp, locus)
        elif item[0] == 'Cellular Component':
            ccdsp = item[1]
            ccdsp = ccdsp.strip()
            ccsym = item[2]
            add2dict2(ccDict, ccsym, ccdsp, locus)
        else:
            print >>sys.stderr, 'Wrong go'
            print >>sys.stderr, go
            print >>sys.stderr, go.split(',')
            sys.exit(0)

#------end addgo------------------------------
def add2dict2(adict, key1, key2, locus):
    if key1 not in adict:
        adict[key1] = set()
    adict[key1].add(locus)
    if key2 not in adict:
        adict[key2] = set()
    adict[key2].add(locus)
################################################
def add2dict(adict, key1, key2, locus):
    if key1 not in adict:
        adict[key1] = []
    adict[key1].append(locus)
    if key2 not in adict:
        adict[key2] = []
    adict[key2].append(locus)
################################################


def readInterpro(interprofile, locusL):
    '''
    pfamDict = {'PF01209':['AT5G57290.3',], 'Ubie':[AT5G57290.3]}
    iprDict = {'IPR001813':['AT5G57290.3',], 'Ribosomal':[AT5G57290.3]}
    mfDict = {'GO:0008168':[AT5G57300.2], 'methy':[AT5G57300.2]}
    bpDict = {'GO:0008168':[AT5G57300.2], 'methy':[AT5G57300.2]}
    ccDict = {'GO:0008168':[AT5G57300.2], 'methy':[AT5G57300.2]}
    '''
    pfamDict = {}
    iprDict = {}
    mfDict = {}
    bpDict = {}
    ccDict = {}
    #pfamsymkeySet = set()
    #pfamdspkeySet = set()
    #iprsymkeySet = set()
    #iprdspkeySet = set()
    for line in open(interprofile):
        lineL = line.rstrip().split('\t')
        locus = lineL[0]
        if locus not in locusL:
            continue
        pfamsym = lineL[4]
        pfamdsp = lineL[5]
        #if 1:
        #    print '%s\t%s' % (pfamsym, pfamdsp)
        #pfamsymkeySet.add(pfamsym)
        #pfamdspkeySet.add(pfamdsp)
        add2dict(pfamDict, pfamsym, pfamdsp, locus)
        iprsym = lineL[11] 
        iprdsp = lineL[12]
        #iprsymkeySet.add(iprsym)
        #iprdspkeySet.add(iprdsp)
        add2dict(iprDict, iprsym, iprdsp, locus)
        if len(lineL) == 14:
            addgo(lineL[13], bpDict, ccDict, mfDict, locus)
        #---End go extract all----------------------------
    #----end read file
    #assert len(pfamsymkeySet) == len(pfamdspkeySet)
    #assert len(iprsymkeySet) == len(iprdspkeySet)
    if TESTdict:
        ct_rdict(pfamDict)
        print >>sys.stderr, "Test finished"
        sys.exit(1)
        
    return (pfamDict, iprDict, bpDict, ccDict, mfDict)
#-------End readInterpro-------------------------------

def output(aDict, keyL, pattern2, locusL, file):
    fl = file + 'Symbol'
    fh = open(fl, 'w')
    print >>fh, "Parameter\t%s" % '\t'.join(locusL)
    count1 = 0
    for item in keyL:
        if pattern2.match(item):
            newLine = [item]
            valueL = aDict[item]
            if isinstance(valueL, list):
                for locus in locusL:
                    newLine.append(str(valueL.count(locus)))  
            elif isinstance(valueL, set):
                for locus in locusL:
                    if locus in valueL:
                        newLine.append('1')
                    else:
                        newLine.append('0')
            print >>fh,  '\t'.join(newLine)
            count1 += 1
        #-----------------
    fh.close()
    #---------finish output pfamsym--------
    fl = file + 'Descrip'
    fh = open(fl, 'w')
    print >>fh, "Parameter\t%s" % '\t'.join(locusL)
    count2 = 0
    for item in keyL:
        if not pattern2.match(item):
            newLine = [item]
            valueL = aDict[item]
            if isinstance(valueL, list):
                for locus in locusL:
                    newLine.append(str(valueL.count(locus)))  
            elif isinstance(valueL, set):
                for locus in locusL:
                    if locus in valueL:
                        newLine.append('1')
                    else:
                        newLine.append('0')
            print >>fh, '\t'.join(newLine)
            count2 += 1
    fh.close()
    #---------finish output pfamsym--------
    return (count1, count2)

def analyze(pfamDict, iprDict, bpDict, ccDict, mfDict, locusL, out):
    pfamKeyL = pfamDict.keys()
    pfamKeyL.sort()
    file = out+'.Pfam'
    (cnt1, cnt2) = \
        output(pfamDict, pfamKeyL, re.compile(r'PF\d{3,}'), locusL, file)
    print >>sys.stderr, "pfamsymbol:%d\npfamdesp:%d" % (cnt1, cnt2)
    iprKeyL = iprDict.keys()
    iprKeyL.sort()
    file= out+'.Interpro'
    (cnt3, cnt4) = \
        output(iprDict, iprKeyL, re.compile(r"IPR\d{4,}"), locusL,
                file)
    print >>sys.stderr, "iprsymbol:%d\niprdesp:%d" % (cnt3, cnt4)
    bpKeyL = bpDict.keys()
    bpKeyL.sort()
    file = out + '.GoBp'
    (cnt5, cnt6) = \
        output(bpDict, bpKeyL, re.compile("GO:"), locusL, file)
    print >>sys.stderr, "bpsymbol:%d\nbpdesp:%d" % (cnt5, cnt6)
    ccKeyL = ccDict.keys()
    ccKeyL.sort()
    file = out + '.GoCc'
    (cnt7, cnt8) = \
        output(ccDict, ccKeyL, re.compile("GO:"), locusL, file)
    print >>sys.stderr, "ccsymbol:%d\nccdesp:%d" % (cnt7, cnt8)
    mfKeyL = mfDict.keys()
    mfKeyL.sort()
    file = out + '.GoMf'
    (cnt9, cnt10) = \
        output(mfDict, mfKeyL, re.compile("GO:"), locusL, file)
    print >>sys.stderr, "mfsymbol:%d\nmfdesp:%d" % (cnt9, cnt10)
    '''
    print >>sys.stderr, 'pfamsym: %d,%d' % (1, cnt1)
    print >>sys.stderr, 'pfamdsp: %d,%d' % (cnt1+1, cnt2+cnt1)
    print >>sys.stderr, 'iprsym:  %d,%d' % (cnt2+cnt1+1,
        cnt2+cnt1+cnt3)
    print >>sys.stderr, 'iprdsp: %d,%d' % (cnt1+cnt2+cnt3+1,\
        cnt1+cnt2+cnt3+cnt4)
    print >>sys.stderr, 'gobpsym:  %d,%d' % (cnt1+cnt2+cnt3+cnt4+1,
        cnt2+cnt1+cnt3+cnt4+cnt5)
    print >>sys.stderr, 'gobpdsp: %d,%d' % (cnt1+cnt2+cnt3+cnt4+cnt5+1,\
        cnt1+cnt2+cnt3+cnt4+cnt5+cnt6)
    print >>sys.stderr, 'goccsym: %d,%d' % \
        (cnt1+cnt2+cnt3+cnt4+cnt5+cnt6+1, \
        cnt1+cnt2+cnt3+cnt4+cnt5+cnt6+cnt7)
    print >>sys.stderr, 'goccdsp: %d,%d' % \
        (cnt1+cnt2+cnt3+cnt4+cnt5+cnt6+cnt7+1, \
        cnt1+cnt2+cnt3+cnt4+cnt5+cnt6+cnt7+cnt8)
    print >>sys.stderr, 'gomfsym: %d,%d' % \
        (cnt1+cnt2+cnt3+cnt4+cnt5+cnt6+cnt7+cnt8+1, \
        cnt1+cnt2+cnt3+cnt4+cnt5+cnt6+cnt7+cnt8+cnt9)
    print >>sys.stderr, 'gomfdsp: %d,%d' % \
        (cnt1+cnt2+cnt3+cnt4+cnt5+cnt6+cnt7+cnt8+cnt9+1, \
        cnt1+cnt2+cnt3+cnt4+cnt5+cnt6+cnt7+cnt8+cnt9+cnt10)
    '''

#----------End analyzeInterpro--------------------
def main():
    print >>sys.stderr, '''
This is used to collect the data for heatmap analysis of proteins
based on domain or go information.    
'''
    print >>sys.stderr, "Print the result to 6 files with the last\
parameter as the prefix "
    if len(sys.argv) != 4:
        print >>sys.stderr, 'Using python %s locus interpro output'\
                % sys.argv[0]
        sys.exit(0)
    #---------------------------------------
    locusL = [line.rstrip() for line in open(sys.argv[1])]
    #if TESTdict:
    #    ct_rlist(locusL)
    #    sys.exit(1)
    (pfamDict, iprDict, bpDict, ccDict, mfDict) = \
            readInterpro(sys.argv[2], locusL)
    analyze(pfamDict, iprDict, bpDict, ccDict, mfDict, locusL,
            sys.argv[3])
if __name__ == '__main__':
    main()

"""
def readInterpro(interprofile):
    '''
    interDict = {
        locus: [[(pfamsym,pfamdsp), ()],
                [(iprsym,iprdsp), ()],
                [set((bpsym, bpdsp),()),
                 set((ccsym, ccdsp),()),
                 set((mfsym, mfdsp),()),
                ]
             ]
    }
    '''
    interDict = {}
    for line in open(interprofile):
        lineL = line.split()
        locus = lineL[0]
        if locus not in interDict:
            interDict[locus] = []
            interDict[locus].append([])
            interDict[locus].append([])
            interDict[locus].append([])
            interDict[locus][2].append(set())
            interDict[locus][2].append(set())
            interDict[locus][2].append(set())
        pfamsym = lineL[4]
        pfamdsp = lineL[5]
        interDict[locus][0].append((pfamsym, pfamdsp))
        iprsym = lineL[11] 
        iprdsp = lineL[12]
        interDict[locus][1].append((iprsym, iprdsp))
        if len(lineL) == 14:
            go = lineL[13]
            for item in go.split(','):
                codonpos = item.find(':')
                leftp = item.find('(')
                rightp = item.find(')')
                if item.startswith('Molecular Function'):
                    mfdsp = item[codonpos:leftp]
                    mfsym = item[leftp+1, rightp]
                    interDict[locus][2][2].add((mfsym, mfdsp))
                elif item.startswith('Biological Process'):
                    bpdsp = item[codonpos:leftp]
                    bpsym = item[leftp+1, rightp]
                    interDict[locus][2][0].add((bpsym, bpdsp))
                elif item.startswith('Cellular Component'):
                    ccdsp = item[codonpos:leftp]
                    ccsym = item[leftp+1, rightp]
                    interDict[locus][2][1].add((ccsym, ccdsp))
                else:
                    print >>sys.stderr, 'Wrong go'
                    sys.exit(0)
            #---End go extract from file-------------------------
        #---End go extract all----------------------------
    #----end read file
    return interDict
#-------End readInterpro-------------------------------

def analyzeInterpro(interDict):
    '''
    interDict = {
        locus: [[(pfamsym,pfamdsp), ()],
                [(iprsym,iprdsp), ()],
                [set((bpsym, bpdsp),()),
                 set((ccsym, ccdsp),()),
                 set((mfsym, mfdsp),()),
                ]
             ]
    }
    '''
"""    
