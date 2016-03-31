#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
#import sys
#sys.path.append("/home/ct/pylib")
'''
This program is used to get the number of how many proteins have repetitons.

Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
==============================================================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=============================================================================================
import re
import sys

if len(sys.argv) != 2:
    print 'Using python %s filename' % sys.argv[0]
    sys.exit(0)


bigger = re.compile('>')
jing = re.compile('#')
try:
    fh = open(sys.argv[1])
except IOError, e:
    print 'can not open file'
else:
    key = 0
    num = 0
    for i in fh:
        if key == 0:
            if bigger.match(i):
                key = 1
        elif key == 1:
            if jing.search(i):
                key = 0
                num += 1
    print 'The number of protein which has repetitions is %d ' % num
    fh.close()
