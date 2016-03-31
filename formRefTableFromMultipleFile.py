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
    print >>sys.stderr, "This is used to construct a reference table \
to show the existence status of terms in each file. The file given \
can be single column file, 1 or 0 would be in output to represent the \
existence of specific term. They can also be two columns file \
separated by tab, the second column is used in output. If one terms \
does not appear in one file, NULL will be used."
    print >>sys.stderr, "Print the result to screen"
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, 'Using python %s multiple files' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    aDict = {}
    fileL = sys.argv[1:]
    termL = []
    none = '0'
    for file in fileL:
        aDict[file] = {}
        for line in open(file):
            lineL = line.strip().split('\t')
            key = lineL[0]
            termL.append(key)
            lenLineL = len(lineL)
            if lenLineL == 1:
                value = '1' 
                none = '0'
            elif lenLineL == 2:
                value = lineL[1]
                none = 'NULL'
            else:
                print >>sys.stderr, "Wrong format"
                sys.exit(1)
            if key not in aDict[file]:
                aDict[file][key]=value
            else:
                print >>sys.stderr, "Duplicate terms in ", file
        #------------Finish reading one file -----------------
    #----------------Finish reading all file--------------
    termL.sort()
    termL = set(termL)
    print "Sample\t%s" % '\t'.join(fileL)
    for term in termL:
        tmpL = [term]
        for file in fileL:
            if term in aDict[file].keys():
                tmpL.append(aDict[file][term])
            else:
                tmpL.append(none)
        #------------------------------------
        print '\t'.join(tmpL)
    #---------------------------------------------------------
#--------------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


