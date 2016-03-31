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

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 5:
        print >>sys.stderr, 'Using python %s allID ctID ykID xjID' % sys.argv[0]
        sys.exit(0)
    '''
    no = 0
    adict = {}
    for file in sys.argv[1:]:
        adict[no] = [line.strip() for line in open(file)]
        no += 1
    #------------------------------------------
    keyL = adict.keys()
    keyL.sort()
    print 'Locus\tDr Wang\tYingKe\tCT'
    for key in keyL:
        for listI in adict[key]:
    '''
    allID = [line.strip() for line in open(sys.argv[1])]
    allID.sort()
    ctID = [line.strip() for line in open(sys.argv[2])]
    ykID = [line.strip() for line in open(sys.argv[3])]
    xjID = [line.strip() for line in open(sys.argv[4])]
    
    print 'Locus\tDr Wang\tYingKe\tCT'
    for id in allID:
        ctflag = 'N'
        ykflag = 'N'
        xjflag = 'N'
        if id in ctID:
            ctflag = 'Y'
            ctID.remove(id)
        if id in ykID:
            ykflag = 'Y'
            ykID.remove(id)
        if id in xjID:
            xjflag = 'Y'
            xjID.remove(id)
        print '%s\t%s\t%s\t%s' % (id, xjflag, ykflag, ctflag)
    #-----------------------------------------

if __name__ == '__main__':
    main()

