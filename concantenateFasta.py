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
timeformat = "%Y-%m-%d %H:%M:%S"
def main():
    introduction='''**This connects each fasta sequence together
with_a given separator. And substitue anycharacter you do not like. 
Print the result to screen.
'''
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, introduction
        print >>sys.stderr, 'Using python %s filename \
delimiter[default X] substitutation["acgtN", "acgtN|ACGT ", default \
uppercaseN] name_for_final_seq[default final]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file  = sys.argv[1]
    if lensysargv > 2:
        delimiter = sys.argv[2]
    else:
        delimiter = 'X'
    if lensysargv > 3:
        formula = sys.argv[3]
    else:
        formula = 'uppercaseN'
    if lensysargv > 4:
        name = sys.argv[4]
    else:
        name = 'final'
    #------------------------------------------
    aDict = {}
    for line in open(file):
        if line[0] == '>':
            locus = line[1:-1]
        else:
            aDict[locus] = line.strip()
    #-------------------------------------------
    seq = delimiter.join(aDict.values())
    if formula == 'uppercaseN':
        seq = seq.upper().replace('N','')
    elif formula.find('|') == -1:
        for i in list(formula):
            seq = seq.replace(i, '')
    elif formula.find('|') != -1:
        formulaL = formula.split('|')
        lena = len(formulaL[0])
        for i in range(lena):
            seq = seq.replace(formulaL[0][i], formulaL[1][i])
    #-----------------------------------------------
    #----------output------------------------
    print name
    print len(seq)
    print seq
#--------------End of main--------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


