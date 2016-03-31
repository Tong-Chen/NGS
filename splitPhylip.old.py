#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, "Split phylip file into given size slide"
        print >>sys.stderr, 'Using python %s phylip size[default 5] \
' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    if lensysargv > 2:
        size = int(sys.argv[2])
    else:
        size = 5
    #----------------------------------------
    phylip = sys.argv[1]
    fh = open(phylip)
    #-------read headline------------
    line = fh.readline()
    species, length = line.split()
        #----form new line for output
    lengthstr = length
    length = int(length)
        #------here floor division needed
    cntSplit = length / size - 1
    last = size + length % size 
    pos = line.rfind(lengthstr)
    assert pos != -1
    blankN = len(lengthstr) - len(str(size))
    newhead = line[:pos] + ' ' * blankN + str(size) + line[pos+len(lengthstr)]
    blankN = len(lengthstr) - len(str(last))
    finalhead = line[:pos] + ' ' * blankN + str(last) + line[pos+len(lengthstr)]
    
#    if species.find(lengthstr) == -1:
#        newhead = line.replace(lengthstr, str(size))
#        finalhead = line.replace(lengthstr, str(last))
#    else:
#        pos = line.rfind(lengthstr)
#        assert pos != -1
#        newhead = line[:pos] + str(size) + line[pos+len(lengthstr)]
#        finalhead = line[:pos] + str(last) + line[pos+len(lengthstr)]
        #---------------------------------------
    species = int(species)
    #-------read first <species> lines
    aDict = {}
    keyL = []
    for i in range(species):
        line = fh.readline()
        lineL = line.split()
        key = lineL[0]
        seq = ' '.join(lineL[1:])
        seq_pos = line.find(seq)
        assert seq_pos != -1
        key = key + ' ' * (seq_pos-len(key)-1)
        aDict[key] = [seq]
        keyL.append(key)
    #-------read following lines---
    while 1:
        #Estimate the end of file and skip blank lines.
        if not fh.readline(): break
        for i in range(species):
            line = fh.readline()
            aDict[keyL[i]].append(line.strip())
    #------------------------------------------------------
    fh.close()
    #-------join sequence---------------
    for key, value in aDict.items():
        aDict[key] = ''.join(value).replace(' ','')
        #instead of replacing blanks in each sequence, here replacing
        #them all together.
    #-------Extract sequence-------------------
    for i in range(cntSplit):
        start = i * size
        output = '.'.join((phylip, str(start)))
        fh = open(output, 'w')
        print >>fh, newhead,
        for item in keyL:
            print >>fh, \
                "%s %s" % (item, aDict[item][start:start+size])
        fh.close()
    #-------Extract final sequence----------------
    start = cntSplit * size
    output = '.'.join((phylip, str(start)))
    fh = open(output, 'w')
    print >>fh, finalhead,
    for item in keyL:
        print >>fh, \
            "%s %s" % (item, aDict[item][start:])
    fh.close()
    #----------------End------------------------------
#-------------End main----------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


