#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
import sys
#sys.path.append("/home/ct/pylib")
'''
This file is used to sort accorderind to one or several columns.
The using is 'ordered.py filename head column[...]'.


Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
==============================================================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=============================================================================================
if len(sys.argv) < 5:
    print 'Using %s filename head spliter column[...]' % sys.argv[0]
    print "%s filename 1 ' ' 2[...]" % sys.argv[0]
    sys.exit(0)

head = sys.argv[2]
spliter = sys.argv[3]
#print '8%s8' % spliter
keycol = int(sys.argv[4]) - 1
if head == '1':
    flag = 0 
adict = {}
#----------------get the sequence and save in a dict----------------------
for line in open(sys.argv[1]):
    if head == '1' and flag == 0:
        flag += 1
        print line,
        continue
    if (head == '1' and flag != 0) or head == '0':
        key = line.split(spliter, 2)[keycol]
        if key in adict:
            adict[key].append(line)
        else:
            adict[key] = [line]
        
#-------------------do the sort and output--------------------------
keylist = adict.keys()
keylist.sort()
#print keylist
for key in keylist:
    for line in adict[key]:
        print line,

if __name__ == '__main__':
    pass

