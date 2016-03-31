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

def main():
    print >>sys.stderr, "\nPrint the result to file, find special\
group contains low complexicity data.\n"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s groupfile complexicity\n' \
            % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------------
    
    lineDict = {}
    savedKey = []
    obseletekeyl = []
    i = -1
    complexicity = float(sys.argv[2])

    for line in open(sys.argv[1]):
        if line[0] == '>':
            savedKey.append([])
            obseletekeyl.append([])
            i += 1
            continue
        elif line[0].isdigit() and line.find('|') != -1:
            key, nouse, seq = line.split('\t')[0:3]
            savedKey[i].append(key)
            lineDict[key] = ''
            
            posDict = {}
            for i40 in seq:
                if posDict.has_key(i40):
                    posDict[i40] += 1
                else:
                    posDict[i40] = 1
            #---------------------------------
            values = posDict.values()
            values.sort()
            #print values
            print values[-1] / len(seq)
            if values[-1] / len(seq) > complexicity:
                obseletekeyl[i].append(key)
        lineDict[key] += line
    #--------End reading-----------------------------------------
    #for i47 in range(i+1):
    #    obseletekeyl[i47] = set(obseletekeyl[i47])
    #------------------------------------------------------------
    #print savedKey
    #print obseletekeyl
    #------------------------------------------------------------
    remaining = 0
    outputfile = sys.argv[1]+'.'+str(complexicity)
    fh = open(outputfile, 'w')
    for i47 in range(i+1):
        savedGroup = savedKey[i47]
        obselGroup = obseletekeyl[i47]
        lenog = len(obselGroup)
        if lenog == 0 or len(savedGroup) - lenog:
            remaining += 1
            print >>fh, '>'
            for keyid in savedGroup:
                if lenog:
                    if keyid not in obselGroup:
                        print >>fh, lineDict[keyid],
                else:
                    print >>fh, lineDict[keyid],
                #------------------------------
            #--End for keyid---------------------
        #---End if lenog-------------------------
    #-------End for i47--------------------------
    fh.close()
    print >>sys.stderr, 'Remaing %d group' % remaining
    #-----------------------------------------------------------    
if __name__ == '__main__':
    main()

