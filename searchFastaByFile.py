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
import re

timeformat = "%Y-%m-%d %H:%M:%S"

def main():
    lensysargv = len(sys.argv)
    if lensysargv < 3:
        print >>sys.stderr, "This will find the position of given \
sequences in a file in FATSA format. The output is a sorted bed file \
Default, the reverse complementary sequence will also be indexed \
case unmatterly. Print the result to files"
        print >>sys.stderr, 'Using python %s fasta seq_file[fasta-format] \
case_matter[default 0, do not consider case. If 1 consider case] \
reverse_complemient[default TRUE, accept 0 means FALSE]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file = sys.argv[1]
    seq_file = sys.argv[2]
    if lensysargv > 3:
        case = sys.argv[3]
    else:
        case = '0'
    if lensysargv > 4:
        rev_com = sys.argv[4]
    else:
        rev_com = 1
    #-----------------------------------------------------------
    patternD = {}
    for i in open(seq_file):
        if i[0] == '>':
            key = i[1:-1]
        else:
            patternD[key] = re.compile(i.strip())
    #--------------------------------------
    lineL = []
    key = ''
    for line in open(file):
        if line[0] == '>':
            if lineL and key:
                seq = ''.join(lineL)
                if not case:
                    seq = seq.upper()
                for k,p in patternD.items():
                    iterator = p.finditer(seq)
                    for iter in iterator:
                        mat = iter.group()
                        start = iter.start()
                        end = start + len(mat)
                        print '\t'.join([key[1:-1],str(start),str(end),mat,
                            k])
                #-------------------------------
            #------------------------------------------- 
            key = line
            lineL = []
        else:
            lineL.append(line.strip())
    #---------------------------------------------------------
    if lineL:
        seq = ''.join(lineL)
        if not case:
            seq = seq.upper()
        for k, p in patternD.items():
            iterator = p.finditer(seq)
            for iter in iterator:
                mat = iter.group()
                start = iter.start()
                end = start + len(mat)
                print '\t'.join([key[1:-1],str(start),str(end),mat, k])
        #-------------------------------
#------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


