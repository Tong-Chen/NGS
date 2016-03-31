#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
import sys
#sys.path.append("/home/ct/pylib")
'''
This is used to get the reverse complement sequences.

Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
==============================================================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=============================================================================================
if len(sys.argv) != 2:
    print 'Using python %s string' % sys.argv[0]
    sys.exit(0)

from string import maketrans, translate

myTranslate = maketrans('AGCTagct', 'TCGAtcga')

print translate(sys.argv[1], myTranslate) 
if __name__ == '__main__':
    pass

