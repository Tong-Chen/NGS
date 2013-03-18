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
    print >>sys.stderr, "Transfer gene annotation bed downloaded from\
ucsc to normal 6-columns bed."
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s geneBed genomeSize' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    aDict = dict([line.strip().split() for line in open(sys.argv[2])])
    alist = [1000,2000,3000,5000]
    #-----------------------------------
    lineno = 0
    for line in open(sys.argv[1]):
        lineno = lineno + 1
        lineL = line.strip().split()
        #print lineL
        chrSta = int(lineL[1])
        chrEnd = int(lineL[2])
        gene = lineL[3] + '_' + str(lineno)
        lineL[3] = gene
        #--------whole gene---------------
        print "\t".join(lineL[:6])
        genesize = chrEnd - chrSta
        #--------exon---------------------
        exonNum   = int(lineL[-3])
        #print '---',lineL[-2],'---'
        #print '---',lineL[-2].split(','),'---'
        exonSize  = [int(i) for i in lineL[-2].split(',') if i]
        exonStart = [int(i) for i in lineL[-1].split(',') if i]
        for i in range(exonNum):
            start = chrSta + exonStart[i] #0-based
            end = start + exonSize[i] #not include
            print "%s\t%d\t%d\t%s\t0\t%s" % (lineL[0], start, end, \
                gene+r'.E', lineL[5])
        #------intron---------------------
        start = 0
        for i in range(exonNum):
            if start < exonStart[i]:
                start = chrSta + start
                end   = chrSta + exonStart[i] 
                print "%s\t%d\t%d\t%s\t0\t%s" % (lineL[0], start, end, \
                    gene+r'.I', lineL[5])
            start = exonStart[i] + exonSize[i]
        if start < genesize:
                start = chrSta + start
                end   = chrSta + genesize
                print "%s\t%d\t%d\t%s\t0\t%s" % (lineL[0], start, end, \
                    gene+r'.I', lineL[5])
        #----------------TSS-----------------------
        theend = int(aDict[lineL[0]])
        #----------------1kb,2kb,5kb---------------
        for item in alist:
            start = 0 if chrSta-item < 0 else chrSta-item
            end   = theend if chrSta+item > theend else chrSta+item
            label = gene + '.' + str(item) + r'.TSS'
            print "%s\t%d\t%d\t%s\t0\t%s" % (lineL[0], start, end,
                    label, lineL[5])
        #--------------TTS-----------------------
        for item in alist:
            start = 0 if chrEnd-item < 0 else chrEnd-item
            end   = theend if chrEnd+item > theend else chrEnd+item
            label = gene + '.' + str(item) + r'.TTS'
            print "%s\t%d\t%d\t%s\t0\t%s" % (lineL[0], start, end,
                    label, lineL[5])
        #----------------------------------------
#---------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


