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



def parse(aDict, locus):
    lena = len(aDict[locus])
    if lena == 1:
        for lineL in aDict[locus][0]:
            lineL[4] = 'm'
            print '\t'.join(lineL)
    else:
        for i in range(lena-1):
            alineL = aDict[locus][i]
            lenAlineL = len(alineL)
            for j in range(i+1, lena):
                blineL = aDict[locus][j]
                lenBlineL = len(blineL)
                for k in range(lenAlineL):
                    epos = alineL[k][3].rfind(r'.E.')
                    if epos > 0:
                        pos = alineL[k][0:3]
                        aDict[locus][i][k][4] = \
                            int(aDict[locus][i][k][4])
                        for l in range(lenBlineL):
                            nepos = blineL[l][3].rfind(r'.E.')
                            if nepos > 0:
                                npos = blineL[l][0:3]
                                aDict[locus][j][l][4] = \
                                  int(aDict[locus][j][l][4])
                                if pos == npos:
                                    aDict[locus][i][k][4] += 1
                                    aDict[locus][j][l][4] += 1
                                    break
                            #-----------------------------
                        #------------------------------------------
                    #-----------------------------------------------
                    else:
                        ipos = alineL[k][3].rfind(r'.I.')
                        if ipos > 0:
                            pos = alineL[k][0:3]
                            for l in range(lenBlineL):
                                nipos = blineL[l][3].rfind(r'.I.')
                                aDict[locus][i][k][4] = \
                                  int(aDict[locus][i][k][4])
                                if nipos > 0:
                                    npos = blineL[l][0:3] 
                                    aDict[locus][j][l][4] = \
                                      int(aDict[locus][j][l][4])
                                    if pos == npos:
                                        aDict[locus][i][k][4] += 1
                                        aDict[locus][j][l][4] += 1
                                        break
                                    #----------------------------
                                #---------------------------------
                            #----------------------------------------
                        #-----------------------------------------
                    #--------------------------------------------
                #----------------------------------------------------
            #--------------------------------------------------------
        #-begion output--------------------------------------------------
        for i in range(lena):
            for lineL in aDict[locus][i]:
                if isinstance(lineL[4], basestring):
                    lineL[4] = 'p'
                else:
                    lineL[4] = str((lineL[4] + 1)/ float(lena))
                print '\t'.join(lineL)
        #------------------------------------------
    #-------------------------------------
#-------------------------------------

#------------------------------------------
def main():
    lensysargv = len(sys.argv)
    if lensysargv != 2:
        print >>sys.stderr, "Detect the exon type and transcript \
type and save in the fifth column. m means no AS. p means with AS.\
Number 1 means constitutive exon or intron. \
Number less than 1 means splicing exon or intron."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s mm9.flat.feature' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    aDict = {}
    endOneG = r'.5000.TTS'
    fh = open(sys.argv[1])
    oldlocus = ''
    while 1:
        line = fh.readline()
        #---end of file--------------
        if not line:
            break
        #---------------------------
        lineL = line.split()
        locus = lineL[3].split('_')[0]
        if oldlocus and locus != oldlocus:
            parse(aDict, oldlocus)
            aDict = {}
        if not aDict:
            i = -1
            aDict[locus] = []
            oldlocus = locus
        #---------------------------------
        i += 1
        aDict[locus].append([])
        aDict[locus][i].append(lineL[:]) 
        #-----------------------------------------
        line = fh.readline()
        while line.find(endOneG) == -1:
            aDict[locus][i].append(line.split()[:])
            line = fh.readline()
        #-----------add last line
        aDict[locus][i].append(line.split()[:])
        #-----------------------------------------------------
    #--------------------------------------
    if aDict:
        parse(aDict, locus)
    fh.close()
#-------------------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


