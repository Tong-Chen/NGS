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
    print >>sys.stderr, "Find common between two files(one \
with_isomer, the other not) and output isomer to the screen"
    if len(sys.argv) != 3:
        print 'Using python %s noiso iso' % sys.argv[0]
        sys.exit(0)
    alist1 = [line.strip() for line in open(sys.argv[1])]
    sys.stdout.writelines([line for line in open(sys.argv[2]) \
            if line.strip().split('.')[0] in alist1])

if __name__ == '__main__':
    main()

