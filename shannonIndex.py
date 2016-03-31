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
from ctIO import readRep
from ctAssist import shannonIndex as si

def main():
    print >>sys.stderr, "Using the average shannonIndex value \
of a group sequences to represent the last entropy."
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    repDict = {}
    readRep(sys.argv[1], repDict)
    
    for valueL in repDict.values():
        for itemD in valueL:
            entropy = 0
            i_valueS = set(itemD.values())
            for item in i_valueS:
                entropy += si(item)
            entropy = entropy / len(i_valueS)
            print "%.2f" % entropy
if __name__ == '__main__':
    main()

