#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''


Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
import sys
from os import system

def extractLen(filename, lenDict):
    lineno = 0
    for line in open(filename):
        if lineno % 3 == 0:
            locus = line[:-1]
        elif lineno % 3 == 1:
            lenDict[locus] = int(line[:-1])
        lineno += 1
#----------End of extract length of protein-------------------------

def extractPos(filename, posDict):
    for line in open(filename):
        if line[0] == '>':
            locus = line[1:-1]
            posDict[locus] = []
        else:
            line = line.rstrip('#\n\r')
            for seqpos in line.split('#'):
                seq, pos = seqpos.split(':')
                posDict[locus].append([int(pos), int(pos)+len(seq)])
        #-------------------------------------------------------------
#-----------End of extractPos-----------------------------------------
def compare(lenDict, posDict, distriDict):
    '''
    posDict = {locus:[[start,end],[start,end]]}
    distriDict = {locus:[[startpercent, endpercent],]}
    '''
    for key, value in posDict.items():
        distriDict[key] = []
        len = lenDict[key]
        for pos in value:
            distriDict[key].append([pos[0]/len, pos[1]/len])
#------------End of compare------------------------------------------

def staticsProDistri(distriDict, startpercent, endpercent, proPerDict):
    '''
    distriDict = {locus:[[startpercent, endpercent],]}
    perDict = {percent:num}
    Attention: contain endpercent, but not startpercent
    (startpercent, endpercent]
    '''
    key0 = '('+str(startpercent)+','+str(endpercent)+']'
    proPerDict[key0] = 0
    
    for key, value in distriDict.items():
        for pospair in value:
            startPos = pospair[0]
            endPos = pospair[1]
            if (((startPos < startpercent) and (endPos > startpercent)\
                    and (endPos+startPos > startpercent * 2))\
                # * -------*--------
                or \
                ((startPos > startpercent) and (endPos <= endpercent))\
                #---*--------*--------
                or \
                ((endpercent > startPos > startpercent) and \
                    (endPos > endpercent) and \
                    (endPos+startPos <= endpercent * 2))):
                #--------*--- *

                proPerDict[key0] += 1
                #break
#------------End of staticsProDistri---------------------------------

def staticsRepDistri(distriDict, startpercent, endpercent, repPerDict):
    '''
    distriDict = {locus:[[startpercent, endpercent],]}
    perDict = {percent:num}
    Attention: contain endpercent, but not startpercent
    (startpercent, endpercent]
    '''
    #key0 = '('+str(startpercent)+','+str(endpercent)+']'
    repPerDict[endpercent] = 0
    
    for key, value in distriDict.items():
        for pospair in value:
            startPos = pospair[0]
            endPos = pospair[1]
            if (((startPos < startpercent) and (endPos > startpercent)\
                    and (endPos+startPos > startpercent * 2))\
                # * -------*--------
                or \
                ((startPos > startpercent) and (endPos <= endpercent))\
                #---*--------*--------
                or \
                ((endpercent > startPos > startpercent) and \
                    (endPos > endpercent) and \
                    (endPos+startPos <= endpercent * 2))):
                #--------*--- *

                repPerDict[endpercent] += 1
#------------End of staticsRepDistri---------------------------------
def outputPerDict(perDict, fh):
    '''
    perDict = {percent:num}
    '''
    keyList = perDict.keys()
    keyList.sort()
    for key in keyList:
        print >>fh, '%.1f\t%d' % (key, perDict[key])
#---------End of outputPerDict---------------------------------------
def drawPic(picname, perDict):
    keyList = perDict.keys()
    valueList = perDict.values()
    keyList.sort()
    valueList.sort()
    step = keyList[0]
    miny = valueList[0] / 2
    maxpy = valueList[-1] + miny
    halfstep = step / 2
    #keyList.insert(0, 0.0)
    #perDict[0] = perDict[keyList[1]]
    tempfile = '.'.join((picname, 'tempct123'))
    fh = open(tempfile, 'w')
    for key in keyList:
        print >>fh, '%.2f\t%d' % (key-halfstep, perDict[key])
    fh.close()
    
    gnu = '.'.join((picname, 'plt'))
    output = '.'.join((picname, 'ps'))
    fh = open(gnu, 'w')
    print >>fh, 'set grid'
    print >>fh, 'set term postscript enhanced color'
    print >>fh, 'set output \"%s\"' % output
    print >>fh, 'set title \"%s repetitions distribute\"' % picname
    print >>fh, 'set xlabel \"regions percent\"'
    print >>fh, 'set ylabel \"# of repetitions\"'
    print >>fh, 'set nokey'
    print >>fh, 'set yrange [%d:%d]' % (miny, maxpy)
    print >>fh, 'set xrange [0:1]'
    print >>fh, 'set xtics %.2f' % step
    print >>fh, 'plot(\"%s\") with boxes' % tempfile
    print >>fh, 'exit'
    fh.close()
    gnuplot = "gnuplot " + gnu
    ps2pdf = "ps2pdf " + output
    system(gnuplot)
    system(ps2pdf)
    if 1:
        rmtemp = "/bin/rm " + tempfile
        rmps = "/bin/rm " + output
        rmgnu = "/bin/rm " + gnu
        system(rmtemp)
        system(rmps)
        system(rmgnu)

#-------------------------------------------------------------------
def main():
    print >>sys.stderr, "Get the distribute data of repetitions"
    if len(sys.argv) != 3 and len(sys.argv) != 6:
        print >>sys.stderr, 'Using python %s forcfile resultfile \
[start(0.0), end(1.0), step(0.1)]' % sys.argv[0]
        sys.exit(0)

    start = 0.0
    end = 0.96
    step = 0.1
    if len(sys.argv) == 6:
        start = float(sys.argv[3])
        end = float(sys.argv[4])
        step = float(sys.argv[5])
    proDistri = sys.argv[2] + '.proDistri.'+str(start)+'.'+str(end)+'.'+str(step)
    lenDict = {}
    posDict = {}
    #posDict = {locus:[[start,end],[start,end]]}
    distriDict = {}
    extractLen(sys.argv[1], lenDict)
    extractPos(sys.argv[2], posDict)
    compare(lenDict, posDict, distriDict)
    #proPerDict = {}
    repPerDict = {}
    while start < end:
        #staticsProDistri(distriDict, start, start+step, proPerDict)
        staticsRepDistri(distriDict, start, start+step, repPerDict)
        start += step
    #outputPerDict(proPerDict)
    #print '---------------------------------------------'
    fh = open(proDistri, 'w')
    outputPerDict(repPerDict, fh)
    fh.close()
    drawPic(proDistri, repPerDict)

if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()')
    main()

