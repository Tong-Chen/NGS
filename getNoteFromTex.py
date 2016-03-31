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
    add = 0
    note = ''
    for line in open(sys.argv[1]):
        start = line.find(r"\note")
        if start != -1:
            add = 1
            if note:
                print note
            note = line[start+5:-1]
        elif add:
            if line.find('}') != -1:
                note += ' ' + line[:line.find('}')]
                add = 0
            else:
                note += ' ' + line.strip()
    #--------------------------
    if note:
        print note

if __name__ == '__main__':
    main()

