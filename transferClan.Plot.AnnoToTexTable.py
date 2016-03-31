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
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    print r'''\begin{table}{lrrl}\footnotesize
\hline
Clan AC & Percentage(\%) & Clan ID & Clan DE \\ \hline'''
    for line in open(sys.argv[1]):
        line = line.replace('_', r'\_')
        lineL = line.strip().split(' ',3)
        print ' & '.join(lineL), r"\\ \hline"
    print r'\end{table}'
if __name__ == '__main__':
    main()

