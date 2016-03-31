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

def main():
    if len(sys.argv) != 3:
        print >>sys.stderr,'Using python %s Pfam-c interpro' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------------------
    clanDict = {} #clanDict = {clan:[pfamlocus]}
    
    for line in open(sys.argv[1]):
        if line.startswith('#=GF AC'):
            clan = line.strip().split()[2]
            if clan in clanDict:
                print >>sys.stderr, 'Wrong'
                sys.exit(1)
            else:
                clanDict[clan] = set()
        elif line.startswith('#=GF MB'):
            pfam = (line.strip().split()[2])[:-1]
            clanDict[clan].add(pfam)
    #-------------------------------------------------
    #print clanDict
    pfamDict = {}
    for line in open(sys.argv[2]):
        lineL = line.split()
        agi = lineL[0]
        pfam = lineL[4]
        if pfam not in pfamDict:
            pfamDict[pfam] = set()
        pfamDict[pfam].add(agi)
    #-------------------------------------
    #print pfamDict

    #--each clan--------------------------------------------
    clannedPfam = set()
    for clan, pfamL in clanDict.items():
        agiList = set()
        for pfam in pfamL:
            if pfam in pfamDict:
                clannedPfam.add(pfam)
                for agi in pfamDict[pfam]: 
                    agiList.add(agi)
        #-------------------------------------
        if len(agiList):
            fh = open(clan, 'w')
            print >>fh, "\n".join(agiList)
            fh.close()
    #--each clan--------------------------------------------
    #---pfam not clanned----------------------------------------
    if len(clannedPfam) != len(pfamDict):
        agiList = set()
        for pfam, agis in pfamDict.items():
            if pfam not in clannedPfam:
                for agi in agis:
                    agiList.add(agi)
        #-----------------------
        if len(agiList):
            fh = open('CL_noclan', 'w')
            print >>fh, '\n'.join(agiList)
            fh.close()
    #---pfam not clanned----------------------------------------

if __name__ == '__main__':
    main()

