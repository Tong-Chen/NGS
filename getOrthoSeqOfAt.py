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
sys.path.append('/home/CT/server/pybin/')
from ctIO import readStdPro

def readOrtholog(file):
    '''
    orthologDict = {
        atitem : set(
            ortholog1, ortholog2, ...        
        )
    }

    paralogDict = {
        atitem : set(
            paralog1, paralog2, ...
        )
    }
    '''
    orthologDict = {}
    paralogDict = {}
    for line in open(file):
        line = line.rstrip()
        lineL = line.split()
        atList = [item for item in lineL if item[:2] == 'AT']
        otherL = [item for item in lineL if item[:2] != 'AT']
        for atitem in atList:
            #---ortholog---------------------
            if atitem not in orthologDict:
                orthologDict[atitem] = set()
            #---ortholog---------------------
            for otherI in otherL:
                orthologDict[atitem].add(otherI)
            #---ortholog---------------------
            #---paralog-----------------------
            if atitem not in paralogDict:
                paralogDict[atitem] = set()
            #---paralog-----------------------
            for x in atList:
                if x != atitem:
                    paralogDict[atitem].add(x)
            #---paralog-----------------------
        #------------End for-----------------
    return (orthologDict, paralogDict)
#------------------------------------------------
def output(fn, orthPara, pepList, repList):
    '''
    orthPara = {
        atitem : set(
            ortholog1, ortholog2, ...        
        ),
    }

    Output format:
    =atitem
    >atitem
    QAZHSDSJKDSDJKDSDJDKSJDJKSDJK
    >ortholog1
    QAZHSDSJKDSDJKDSDJDKSJDJKSDJK
    >ortholog2
    QAZHSDSJKDSDJKDSDJDKSJDJKSDJK
    '''
    fh = open(fn, 'w')
    old = sys.stdout
    sys.stdout = fh
    orthParaKL = orthPara.keys()
    orthParaKL.sort()
    for orthParaK in orthParaKL:
        #------modified 2011-06-27----only output proteins have
        #repetitions
        if orthParaK not in repList:
            continue
        #------modified 2011-06-27----only output proteins have
        #repetitions
        orthParaV = orthPara[orthParaK]
        num = len(orthParaV) + 1
        #------modified 2011-06-28----delete those proteins
        #with_only one ortholog or paralog which is itself.
        if num == 1:
            continue
        #------modified 2011-06-28----delete those proteins
        #with_only one ortholog or paralog which is itself.
        print '=%s %d' % (orthParaK, num) 
        #---------At sequence----------------------
        print '>%s\n%s' % \
            (orthParaK, (pepList[-1])[orthParaK])
        #--------ortholog sequence----------------
        for orthParaId in orthPara[orthParaK]:
            noThisId = 1
            for aDict in pepList:
                if orthParaId in aDict:
                    print '>%s\n%s' % \
                        (orthParaId, aDict[orthParaId])
                    noThisId = 0
                #---if it belongs to this dict------
            #---------finish trace--------------
            if noThisId:
                print >>sys.stderr, "No %s" % orthParaId
                sys.exit(1)
    #--End one ortholog group----------------------------
        
    fh.close()
    sys.stdout = old
#-----------------------------------------------

def main():
    print >>sys.stderr, "Print the result to file"
    if len(sys.argv) != 5:
        print >>sys.stderr, 'Using python %s orthologs pepPath\
outputPrefix replocus' % sys.argv[0]
        sys.exit(0)
    #----------------------------------------------
    #-------read ortholog paralog------------------
    orthoD, paraD = readOrtholog(sys.argv[1])
    #-------read ortholog paralog------------------
    pepfile = ['Bdistachyon_114_peptide.fa','Cpapaya_113_peptide.fa',
            'osa1r6.1_pep','Ptrichocarpa_156_peptide.fa',
            'Sbicolor_79_peptide.fa', 'Vvinifera_145_peptide.fa',
            'Zmays_121_peptide.fa', 'TAIR_PEP']
    
    pepList = ['BdiD', 'CpaD', 'OsaD', 'PtrD', 'SbiD',
            'VviD', 'ZmaD', 'AtD']
    
    #---------Read Seq-----------------------------------------
    path = sys.argv[2]
    if path[-1] != '/':
        path += '/'
    for i in range(8):
        pathpepfile = path+pepfile[i]
        pepList[i] = readStdPro(pathpepfile)
    #---------Read Seq-----------------------------------------
    repList = [locus.strip() for locus in open(sys.argv[4])]
    #---------Output ortholog----------------------------------
    fn = sys.argv[3] + '.Ortholog.For.C'
    output(fn, orthoD, pepList, repList)
    #---------Output ortholog----------------------------------
    fn = sys.argv[3] + '.Paralog.For.C'
    output(fn, paraD, pepList, repList)


if __name__ == '__main__':
    main()

