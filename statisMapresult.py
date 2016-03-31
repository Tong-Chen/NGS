#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
import sys

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    aDict = {}
    '''
    aDict = {At : [total ortho/paralog, # of map, # of unmap}
    '''
    for line in open(sys.argv[1]):
        if line[0] == '=':
            locus, num = line[1:-1].split()
            num = int(num)
            map = 0
            unmap = 0
            aDict[locus] = [num, map, unmap]
        elif line[0] != '>' and (line.find('#') != -1 or 
            line.find(':') != -1) and \
            line.find('-') == -1:
            (aDict[locus])[1] += 1
        elif line[0] == '-':
            (aDict[locus])[2] += 1
    #-------------------------------------------
    
    #print aDict

    for locus, value in aDict.items():
        assert value[1] + value[2] == value[0]
        if value[0] == 1:
            continue
        #changed 2011-08-16
        #bit = value[1] / value[0] #old operation
        bit = (value[1]-1) / (value[0]-1)
        #changed 2011-08-16
        #print "%.2f" % bit #old operation 
        print "%.2f" % bit 
            
if __name__ == '__main__':
    main()

