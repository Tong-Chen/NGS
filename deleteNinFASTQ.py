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
        print >>sys.stderr, "Delete N in fastq sequences."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename \
begin[default yes] end[default yes] \
inner[default no]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    if lensysargv > 2:
        begin = sys.argv[2]
    else:
        begin = 'yes'
    if lensysargv > 3:
        end = sys.argv[3]
    else:
        end = 'yes'
    if lensysargv > 4:
        inner = sys.argv[4]
    else:
        inner = 'no'
    #--------------------------------------
    delPosL = []
    j = 1
    for line in open(sys.argv[1]):
        if j % 4 == 2:
            delPosL = []
            line = line.strip()
            lenline = len(line)
            #--------delete begin-------------------------
            if begin == 'yes' and inner == 'no':
                for i in range(lenline):
                    if line[i] != 'N':
                        break
                    delPosL.append(i)
            #----------delete end--------------------------
            if end == 'yes' and inner == 'no':
                for i in range(lenline-1,-1,-1):
                    if line[i] != 'N':
                        break
                    delPosL.append(i)
            #---------delete all-----------------------------
            if begin == 'yes' and end == 'yes' and inner == 'yes':
                for i in range(lenline):
                    if line[i] == 'N':
                        delPosL.append(i)
            #--------begin deleting---------------------------------------
            print delByPos(line, delPosL, lenline)
            #print delPosL, lenline
        elif j % 4 == 0:
            line = line.strip()
            print delByPos(line, delPosL, len(line))
        else:
            print line,
        j += 1
    #------------------end reading-----------------------
#---------------------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


