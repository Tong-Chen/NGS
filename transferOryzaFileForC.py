#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
#import sys
#sys.path.append("/home/ct/pylib")
'''
this file is used to transfer oryza for C use.

Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
==============================================================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=============================================================================================
if __name__ == '__main__':
    import sys
    print >>sys.stderr, 'this used to transfer oryza sequence for c\n\
+   and output the result to stdout'

    if len(sys.argv) < 2:
        print >>sys.stderr, 'Using %s filename' % sys.argv[0]
        sys.exit(1)

    flag = 0
    str = ''
    for line in open(sys.argv[1]):
        if line.find('>') == 0:
            key = line[1:7]
            if str != '':
                print len(str)
                print str
                str = ''
            print key
        else:
            str += line[:-1]
    if str != '':
        print len(str)
        print str
        str = ''

