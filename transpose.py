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
    if len(sys.argv) != 2:
        print >>sys.stderr, "Transpos a file, turn rows to columns and\
 columns to rows."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename(- represent STDIN)' % sys.argv[0]
        sys.exit(0)
    #--------------------------------------------------
    file = sys.argv[1]
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    list = []
    [list.append(line.strip().split('\t')) for line in fh]
    column = len(list[0])
    row = len(list)
    for i in range(column):
        tmplist = [list[j][i] for j in range(row)]
        print '\t'.join(tmplist)
    if file != '-':
        fh.close()
#-----------------------------------------------
if __name__ == '__main__':
    main()

