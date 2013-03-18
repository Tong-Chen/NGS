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
###ACTGTAGGCACCATCAATC
import sys
from Levenshtein import ratio,distance

def cliA(line, lenline, adapt, lenadapt, short_mm_len, extend,
        editDratio, shortReadsLen):
    reads = line
    #-----Begin full adaptor-----------------------
    tmpAdapt = adapt
    exactPos = line.rfind(tmpAdapt)
    if exactPos != -1 and exactPos > shortReadsLen:
        if extend == 'yes':
            return line[:exactPos]
        else:
            if exactPos + lenadapt == lenline:
                return line[:exactPos]
        #---------------------------------------------
    else:
        if extend == 'yes':
            end = lenadapt + shortReadsLen - 1 # start, no less than
                                               # shortReadsLen
        else:
            end = lenline -1
        for j in range(lenline, end, -1):
            tmpReads = line[j-lenadapt:j]
            if ratio(tmpReads, tmpAdapt) >= editDratio:
            # and tmpReads[0] == tmpAdapt[0]:
                return line[:j-lenadapt]
        #------------mismatch----------------
    #---------------End full adaptor----------------
    #--------------Begin partly clipped----------
    for i in range(lenadapt-1, short_mm_len, -1):
        tmpAdapt = adapt[:i]
        exactPos = line.rfind(tmpAdapt)
        if exactPos != -1 and exactPos+i == lenline:
            return line[:exactPos]
        else:
            tmpReads = line[lenline-i:lenline]
            if ratio(tmpReads, tmpAdapt) >= editDratio:
                #and \
                #tmpReads[0] == tmpAdapt[0]:
                return line[:lenline-i]
            #------------mismatch----------------
        #---------------End full adaptor----------------
    #--no process-----------------------------------
    return reads
#-----------------------------------------------------------------------

def main():
    len_argv = len(sys.argv)
    if len_argv < 3:
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s fastq adapter \
shortest-match-len[1] extend[[yes].if "yes" delete even adaptor \
located not the end; if "no" delete only when adaptor locates at the end. \
Only suitable for full-length adaptor. ] \
edit_ratio, default 0.85.  similarity no less than given value) \
minimum_length_for_reads_keeping[default 30]' % sys.argv[0]
        sys.exit(0)
    #---------------------------------------------------------------
    fastq = sys.argv[1]
    adapt = sys.argv[2]
    lenadapt = len(adapt)
    if len_argv > 3:
        short_mm_len = int(sys.argv[3]) - 1 #minus one is for program
                                            #usage, no influention to
                                            #result.
    else:
        short_mm_len = 0 #(means at least 1 nucleotide)
    if len_argv > 4:
        extend = sys.argv[4]
    else:
        extend = 'yes'
    if len_argv > 5:
        editDratio = float(sys.argv[5])
    else:
        editDratio = 0.85
    if len_argv > 6:
        shortReadsLen = int(sys.argv[6])
    else:
        shortReadsLen = 30
    #---------------------------------------------
    aDict = {}
    i = 1
    output = 0
    for line in open(fastq):
        line = line.rstrip()
        lenline = len(line)
        if i % 4 == 1:
            lastline = line
        elif i % 4 == 2:
            line = cliA(line, lenline, adapt, lenadapt, short_mm_len, extend,
                    editDratio, shortReadsLen)
            newlen = len(line)
            #------for statistics-------------------
            del_len = lenline - newlen
            if del_len not in aDict:
                aDict[del_len] = 1
            else:
                aDict[del_len] += 1
            #------for statistics-------------------
            #--------------------------------------
            if newlen >= shortReadsLen:
                output = 1
            if output:
                print lastline #the first line name line
                print line  #second,  reads line
        elif i % 4 == 0:
            if output:
                print line[:newlen] #the last line
                output = 0
        #---output the third line---------------------------------
        elif output:
            print line
        i += 1
    #--------------------END of file--------------------------------
    #--------output statistics--------------------------
    fh = open(fastq+".clip.sta", 'w')
    keyL = aDict.keys()
    keyL.sort()
    for key in keyL:
        print >>fh, "%d\t%d" % (key, aDict[key])
    fh.close()
#--------------------------------------------------------------------
if __name__ == '__main__':
    main()

