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

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s groupfile' % sys.argv[0]
        sys.exit(0)
    #---------------------------------------------------------
    fh = open('tmp20110117', 'w')
    for line in open(sys.argv[1]):
        if line[0].isdigit() and line.find('|') != -1:
            list=line.split('\t')
            littlelist = list[3].split(' [')
            num = littlelist[1].count(',')
            if int(littlelist[0]) != num + 1:
                littlelist[0] = str(num + 1)
            list[3] = ' ['.join(littlelist)
            line = '\t'.join(list)
        print >>fh, line,
    fh.close()
    mv = 'mv tmp20110117 '+sys.argv[1]
    os.system(mv)

if __name__ == '__main__':
    main()

