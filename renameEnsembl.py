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
import re

def main():
    print >>sys.stderr, "This used to rename ensemble file."
    
    pattern = re.compile('(.).+?_([a-z]+?)\..+?(pep)\.all\.(.+)')

    filenamelist = os.listdir(os.curdir)
    for file in filenamelist:
        if file.find('_') == -1:
            continue
        match = pattern.match(file)
        if match:
            newname='.'.join((match.group(1), match.group(2), \
                match.group(3), match.group(4)))
            os.renames(file, newname)
        else:
            print >>sys.stderr, file, 'wrong'
            sys.exit(0)
if __name__ == '__main__':
    main()

