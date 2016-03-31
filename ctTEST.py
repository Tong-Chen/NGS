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
from random import choice

def ct_rdict(aDict, cnt=10):
    keyL = aDict.keys()
    for i in range(cnt):
        key = choice(keyL)
        print key, aDict[key]
    
def ct_rlist(aList, cnt=10):
    for i in range(cnt):
        print choice(aList)


