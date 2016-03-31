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

def delByPos(line, delPosL, lenline):
    return ''.join([line[i] for i in range(lenline) if i not in
        delPosL])
#--------------------------------------------------


def main():
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, "Delete PolyA in fastq sequences."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename[- represents \
stdin] string[AAAAA, default at least 5A]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file = sys.argv[1]
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    if lensysargv > 2:
        string = sys.argv[2]
    else:
        string = 'AAAAA'
    #--------------------------------------
    lenstr = len(string)
    delPosL = []
    j = 1
    for line in fh:
        if j % 4 == 2:
            delPosL = []
            line = line.strip()
            lenline = len(line)
            posA = line.rfind(string)
            if posA+lenstr == lenline:
                for lenline in range(posA-1,30,-1):
                    if line[lenline] != 'A' and line[lenline] != 'N':
                        lenline += 1
                        break
            #-----------------------------------
            print line[:lenline]
        elif j % 4 == 0:
            line = line.strip()
            print line[:lenline]
        else:
            print line,
        j += 1
    #------------------end reading-----------------------
    if file != '-':
        fh.close()
#---------------------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


