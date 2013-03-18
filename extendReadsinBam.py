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
import re
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
def main():
    lensysargv = len(sys.argv)
    if lensysargv < 3:
        print >>sys.stderr, "This is used  to extend mapped reads to given \
length. Only for single end reads and flag is decimal. Print the result to screen."
        print >>sys.stderr, 'Using python %s filename [sam file or  represents \
STDIN] expected_length' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    length_op = ['M', "I", "S", "=", "X"]
    file = sys.argv[1]
    exp_len = int(sys.argv[2])
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------begin reading--------------------
    for line in fh:
        if line[0] == "@":
            sys.stdout.write(line) 
            continue
        lineL = line.split('\t', 11)
        seq   = lineL[9]
        len_seq = len(seq)
        diff = exp_len - len_seq
        #print diff
        if diff <= 0:
            sys.stdout.write(line) 
            continue
        #---------------------------------
        flag  = lineL[1]
        cigar = lineL[5]
        quality = lineL[10]
        #---get numbers and alphabets--------------
        cigarL = []
        numL = [] 
        for i in cigar:
            if i.isdigit():
                numL.append(i)
            elif i.isalpha():
                cigarL.append(''.join(numL))
                cigarL.append(i)
                numL = []
        #-----------------------------------
        #------correct cigar for extension------------------
        #print cigarL
        if flag == '0':
            assert cigarL[-1] == 'M', cigar
            cigarL[-2] = str(int(cigarL[-2]) + diff)
        elif flag == '16':
            assert cigarL[1]  == "M", cigar 
            cigarL[0] = str(int(cigarL[0]) + diff)
            lineL[3] = str(int(lineL[3])-diff) 
        #print cigarL
        #---------------------------------------------------
        #---------extend SEQ and QUAL------------
        newseq = seq + seq[:diff]
        while len(newseq) < 100:
            newseq += newseq
        lineL[9]  = newseq[:100]
        newquality = quality + quality[:diff]
        while len(newquality) < 100:
            newquality += newquality
        lineL[10] = newquality[:100]
        lineL[5]  = ''.join(cigarL)
        sys.stdout.write('\t'.join(lineL)) 
    #---------turn off file-------------------
    if file != '-':
        fh.close()
#-------------------------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


