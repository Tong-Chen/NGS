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
    print >>sys.stderr, "\nPrint the result to screen, find special\
group with or without given go term.\n"
    if len(sys.argv) != 3:
        print >>sys.stderr, "grep 'cell wall' arab.rep.group.anno | sort\
-u | grep 'GO:' | cut -f 2 | sed 's/:$//'\n\n"
        print >>sys.stderr, 'Using python %s groupfile gofile\n' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------------
    goidlist = [goid[3:-1] for goid in open(sys.argv[2])]
    
    lineDict = {}
    savedKey = []
    obseletekeyl = []
    i = -1

    for line in open(sys.argv[1]):
        if line[0] == '>':
            savedKey.append([])
            obseletekeyl.append([])
            i += 1
            continue
        elif line[0].isdigit() and line.find('|') != -1:
            key = line.split('\t',1)[0]
            savedKey[i].append(key)
            keyin = 0
            lineDict[key] = ''
        elif (not keyin) and line.find('GO:') != -1 \
                and line.split(':')[-2] in goidlist:
            obseletekeyl[i].append(key)
            keyin = 1
        lineDict[key] += line
    #--------End reading-----------------------------------------
    #for i47 in range(i+1):
    #    obseletekeyl[i47] = set(obseletekeyl[i47])
    #------------------------------------------------------------
    #print savedKey
    #print obseletekeyl
    #------------------------------------------------------------
    remaining = 0
    for i47 in range(i+1):
        savedGroup = savedKey[i47]
        obselGroup = obseletekeyl[i47]
        lenog = len(obselGroup)
        if lenog == 0 or len(savedGroup) - lenog:
            remaining += 1
            print '>'
            for keyid in savedGroup:
                if lenog:
                    if keyid not in obselGroup:
                        print lineDict[keyid],
                else:
                    print lineDict[keyid],
                #------------------------------
            #--End for keyid---------------------
        #---End if lenog-------------------------
    #-------End for i47--------------------------
    print >>sys.stderr, 'Remaing %d group' % remaining
    #-----------------------------------------------------------    
if __name__ == '__main__':
    main()

