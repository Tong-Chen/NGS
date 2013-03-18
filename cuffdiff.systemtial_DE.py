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
Filename:TP16_TP22_TP23.gene_exp.diff_g.expr_1.0.len_300.TP16___TP22 
(This lists genes which are differentially expressed at two samples
which are labeled by the suffix of filename.
#Only first three columns are needed. Other columns will be ignored. 
One Header line is needed. )
test_id TP22    TP23    TP16    TP22___TP23 TP16___TP22 TP16___TP23
1500015O10Rik__XLOC_000123__=   0.146503    0.718229    26.9616 no  yes yes
1700123I01Rik__XLOC_012935__j   4.0776  0.0541536   0.0254863   no  yes no
2510009E07Rik__XLOC_010388__=   0.217675    9.46837 34.0829 yes yes no

'''


if False:
    print "This program does not work under python 3, \
run in python 2.x."




import sys
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"


def diff(fileL, sampL):
    sampS = list(set([j for i in sampL for j in i]))
    prefix = '.'.join(fileL[0].split('.')[:-1])
    #---common than --------------------
    for samp in sampS:
        commonThanL = []
        specialHighL = []
        otherSampL = sampS[:]
        otherSampL.remove(samp)
        otherSamps = '___'.join(otherSampL)
        for file in fileL:
            sampL = file.split('.')[-1].split('___')
            if samp in sampL:
                #------------------special high -------------
                header = 1
                tmpL = []
                for line in open(file+'.'+samp):
                    if header:
                        header -= 1
                        continue
                    tmpL.append(line.split('\t')[0])
                specialHighL.append(set(tmpL))
                #---------------common than -------------------
                sampL.remove(samp)
                anotherSamp = sampL[0]
                header = 1
                tmpL = []
                for line in open(file+'.'+anotherSamp):
                    if header:
                        header -= 1
                        continue
                    tmpL.append(line.split('\t')[0])
                commonThanL.append(set(tmpL))
            #-----------------------------------------------------
        #---------------END iter all files for one sample---------
        newSet = specialHighL[0]
        for i in specialHighL[1:]:
            newSet = newSet.intersection(i)
        fh = open(prefix+'.'+otherSamps+'.allLowerThan'+samp, 'w')
        newSet = list(newSet)
        newSet.sort()
        for i in newSet:
            print >>fh, i
        fh.close()
        #-------------------------------------
        newSet = commonThanL[0]
        for i in commonThanL[1:]:
            newSet = newSet.intersection(i)
        fh = open(prefix+'.'+otherSamps+'.allHigherThan'+samp, 'w')
        newSet = list(newSet)
        newSet.sort()
        for i in newSet:
            print >>fh, i
        fh.close()
    #------------------END all sample---------------
#----------------diff------------------------------

def main():
    lensysargv = len(sys.argv)
    if lensysargv < 4:
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename_\
at_least_three[unlimited, a series of files. Like \
TP16_TP22_TP23.gene_exp.diff_g.expr_1.0.len_300.TP16___TP22. \
#The suffix of filename must be samp1___samp2. ANs samp1 and samp2 must \
#appear at the second and third column if first line in file.]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    fileL = sys.argv[1:]
    sampL = []
    #-----------end close fh-----------
    fh1 = ''
    fh2 = ''
    for file in fileL:
        header = 1
        if fh1 or fh2:
            fh1.close()
            fh2.close()
        for line in open(file):
            if header:
                samp1,samp2 = file.split('.')[-1].split('___')
                lineL = line.split("\t")
                lenline = len(lineL)
                col1 = ''
                col2 = ''
                for i in range(1,lenline):
                    if lineL[i] == samp1:
                        col1 = i
                    elif lineL[i] == samp2:
                        col2 = i
                    if col1 and col2:
                        break
                #-----------------------------
                assert col1 and col2, file
                sampL.append((samp1, samp2))
               # assert (name1==samp1 and name2==samp2) \
               #     or (name1==samp2 and name2==samp1), file
                fh1 = open(file+'.'+samp1, 'w')
                print >>fh1, '\t'.join([lineL[0],lineL[col1],lineL[col2]])
                fh2 = open(file+'.'+samp2, 'w')
                print >>fh2, '\t'.join([lineL[0],lineL[col1],lineL[col2]])
                header = 0
            else:
                lineL = line.split('\t')
                if (float(lineL[col1])>float(lineL[col2])):
                    print >>fh1,'\t'.join([lineL[0],lineL[col1],lineL[col2]])
                else:
                    print >>fh2,'\t'.join([lineL[0],lineL[col1],lineL[col2]])
            #---------------------------------------------
        #---------------END one file--------------
    #---------------END iterating files----- 
    if fh1 or fh2:
        fh1.close()
        fh2.close()
    #--------------close the last file----------
    diff(fileL,sampL)
#-----------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


