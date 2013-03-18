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
import random
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
def main():
    lensysargv = len(sys.argv)
    if lensysargv != 6:
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s chromosomesize \
number_of_randome_pos number_of_replicate head include_random[1,0]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file = sys.argv[1]
    num_random_pos = int(sys.argv[2])
    num_random_rep = int(sys.argv[3])
    include_random = int(sys.argv[4])
    head = int(sys.argv[5])

    aDict = {}
    for line in open(file):
        if head:
            head -= 1
            continue
        key, length = line.split()
        if (not include_random) and (key.find('random') == -1):
            aDict[key] = int(length) - 1
        elif include_random:
            aDict[key] = int(length) - 1
    #------------------------------------------------------------
    keyL = aDict.keys()
    for i in range(num_random_rep):
        file_op = file + '.' + str(i) + '-' + str(num_random_pos) 
        fh = open(file_op, 'w')
        for j in range(num_random_pos):
            rdm_chr = random.choice(keyL)
            rdm_pos = random.randint(0,aDict[rdm_chr])
            print >>fh, "%s\t%d\t%d\t%d" % (rdm_chr, rdm_pos, rdm_pos+1,j)
        fh.close()
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


