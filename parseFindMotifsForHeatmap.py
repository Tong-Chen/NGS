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
    len_argv = len(sys.argv)
    if len_argv < 4:
        print >>sys.stderr, "Print the result to file"
        print >>sys.stderr, 'Using python %s filename min max head[1]\
strand_info[consider stand info if 1 is given else do not consider \
strand info when 0 is given.]' % sys.argv[0]
        sys.exit(0)
    #--------------------------------------------------------------
#'''
#Offset : Count form left to right no matter + or -
#'''
    if len_argv > 4:
        head = int(sys.argv[4])
    else:
        head = 1
    if len_argv > 5:
        strand = sys.argv[5]
    else:
        strand = 0
    #--------------------------------------
    aDict = {}
    min = int(sys.argv[2])
    max = int(sys.argv[3])
    file = sys.argv[1]
    for line in open(file):
        if head:
            head -= 1
            continue
        #-----------------------
        lineL = line.split("\t")
        start = int(lineL[1])
        end = start + len(lineL[2]) #not inciluded
        if (lineL[0] not in aDict):
            aDict[lineL[0]] = [(start, end, lineL[4])]
        else:
            aDict[lineL[0]].append((start, end, lineL[4]))
    #------------------------------------------------------------
    if strand:
        vDict = {'+':'1',  '-':'-1'}
    else:
        vDict = {'+':'1',  '-':'1'}
    step = 1
    print "Peak\t%s" % '\t'.join([str(i) for i in range(min, max+1, step)])
    for key,  valueL in aDict.items():
        newdict = {}
        for itemL in valueL:
            for i in range(itemL[0], itemL[1]):
                newdict[i] = vDict[itemL[2]]
        #----------------------------------------
        tmpvalue = [key]
        for i in range(min, max+1, step):
            if i not in newdict:
                tmpvalue.append('0')
            else:
                tmpvalue.append(newdict[i])
        #-------------------------------------------
        print '\t'.join(tmpvalue)
        

#-----------------------------------------------------------------------
if __name__ == '__main__':
    main()

