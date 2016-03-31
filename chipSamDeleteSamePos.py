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

def readSam(file):
    '''
    sam file is BWA output, major consideration is the first five
    ccolumns.

    aDict = {chr\tpos:[line,]}
    '''
    header = []
    aDict = {} 
    unmapped = []
    for line in open(file):
        if line[0] == '@':
            header.append(line)
        else:
            lineL = line.split('\t',4)
            if lineL[1] == '4':
                try:
                    assert lineL[2] == '*'
                    assert lineL[3] == '0'
                except AssertionError:
                    print >>sys.stderr, "Wrong with unmapped reads,\
                    %s" % line
                    sys.exit(1)
                #-----------------------------
                unmapped.append(line)
                continue
            #---------------------------------
            key = '\t'.join(lineL[2:4])
            if key in aDict:
                aDict[key].append(line)
            else:
                aDict[key] = [line]
            #--------------------------
        #-----------------------------------
    return header, aDict, unmapped
#---------------------------------------------

def deleteDuplicate(aDict):
    for key, valueL in aDict.items():
        len_v = len(valueL)
        if len_v > 1:
            valueL.sort(key=lambda x: int(x.split('\t',5)[4]),\
                reverse=True)
            aDict[key] = valueL[0]
        #-----------------------------------
        elif len_v:
            aDict[key] = valueL[0] 
        else:
            print >>sys.stderr, "Error of %s" % key
#------------------------------------------
def outputUnmapped(header, unmapped, outputUFile):
    fh = open(outputUFile, 'w')
    print >>fh, ''.join(header),
    print >>fh, ''.join(unmapped),
    fh.close()
#-------------------------------------------------
def outputS(header, aDict, outputSFile):
    fh = open(outputSFile, 'w')
    print >>fh, ''.join(header),
    valueL = aDict.values()
    valueL.sort(key=lambda x: \
        int((x.split('\t',1)[0]).split('.')[-1]))
    fh.close()
#------------------------------------------------------

def main():
    print >>sys.stderr, "Print the result to files"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    input = sys.argv[1]
    outputUFile = input[:-3]+'.Unmapped.sam'
    outputSFile = input[:-3]+'.DelSamePos.sam'
    header, aDict, unmapped = readSam(sys.argv[1])
    deleteDuplicate(aDict)
    outputS(header, aDict, outputSFile)
    outputUnmapped(header, unmapped, outputUFile)
#-----------------------------------------
if __name__ == '__main__':
    from time import localtime 
    startTime = localtime()
    main()
    endTime = localtime()
    fh = open('python.log', 'a')
    print >>fh, "Run time : %s - %s " % (startTime, endTime)
    fh.close()

