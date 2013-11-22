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
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
def main():
    print >>sys.stderr, "Print the result to screen"
    lensysargv = len(sys.argv)
    if lensysargv < 3:
        print >>sys.stderr, 'Using python %s function[comm,firstSpecial] \
outputDir/[make this dir by yourself] \
filenames[common suffix]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    aDict = {'comm':'diffCom.py', 'firstSpecial':'firstSpecial.py',\
            '.comm':'diffCom.py', '.firstSpecial':'firstSpecial.py'}
    function_N = sys.argv[1] 
    function   = aDict[function_N]
    outputDir = sys.argv[2]
    if outputDir[-1] != '/':
        outputDir += '/'
    filenames = sys.argv[3:]
    suffix = filenames[0].split('.')[-1]
    lenfilenames = len(filenames)
    for i in range(lenfilenames-1):
        fi = filenames[i]
        fip = fi.rsplit('.',1)[0]
        for j in range(i+1,lenfilenames):
            fj = filenames[j]
            fjp = fj.rsplit('.',1)[0]
            fijp = outputDir+fip+'_'+fjp+'.'+function_N+'.'+suffix
            #if function_N == 'comm':
            cmd = ' '.join((function, fi, fj,'|sort -u ','>',fijp,'&')) 
            print cmd
            os.system(cmd)
            if function_N == 'firstSpecial' or function_N == '.firstSpecial':
                fjip = outputDir+fjp+'_'+fip+'.'+function_N+'.'+suffix
                cmd = ' '.join((function, fj, fi,'|sort -u','>',fjip,'&')) 
                print cmd
                os.system(cmd)
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


