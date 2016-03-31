#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''


Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
#from __future__ import division, with_statement
import sys

def main():
    print >>sys.stderr, "Get locus related group"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s locusfile groupfile' % sys.argv[0]
        sys.exit(0)
    
    #groupFlag = 0
    groupDict = {}
    for line in open(sys.argv[2]):
        if line[0] == '>':
            groupId = line
            groupDict[groupId] = ''
        else:
            groupDict[groupId] += line
    
    for locus in open(sys.argv[1]):
        locus = locus[:-1]
        for groupId, groupInfo in groupDict.items():
            if locus in groupInfo:
                print groupId
                print groupInfo
                del groupDict[groupId]
if __name__ == '__main__':
    main()

