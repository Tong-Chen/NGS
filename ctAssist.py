#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================

def shannonIndex(pep, num=20):
    '''
    shannonIndex += Pi * log2(Pi) for Pi in num of elements.
    Pi = (num of i) / (lenpep)
    '''
    from math import log
    adict = {}
    for i in pep:
        if i in adict:
            adict[i] += 1
        else:
            adict[i] = 1
    #----------------------------------
    entropy = 0
    log2 = log(2)
    lenpep = len(pep)
    for count in adict.values():
        pi = count / lenpep
        entropy += pi * log(pi) / log2
    return (-1) * entropy

#alphabet = [chr(i) for i in range(97,123)]

def labelP(num, alphabet):
    division = num / 26
    remain = num % 26
    label = ''
    label = alphabet[remain] + label
    while division > 25:
        remain = division % 26
        division = division / 26
        label = alphabet[remain-1] + label
    #---------------------------
    if division > 0:
        label = alphabet[division-1] + label
    return label
#-----------------------------------


if __name__ == '__main__':
    print shannonIndex("aaaaaaaaaa")
    print shannonIndex("ABCDEFGHIJKLMNOPQRSTABCDEFGHIJKLMNOPQRST")
    str = 'V' * 16 + 'A' * 2
    print str, shannonIndex(str)
    str = 'L' + 'T' * 17
    print str, shannonIndex(str)
    print shannonIndex("ababababab")

    str = 'M' * 9 + 'L' * 8 + '-'
    print str, shannonIndex(str)

    str = 'A' + 'S'*7 + '-'*10
    print str, shannonIndex(str)
    str = 'A' *3 + '-' + 'L'*14
    print str, shannonIndex(str)
    alphabet = [chr(i) for i in range(97,123)]
    for i in range(52):
        print labelP(i, alphabet)
        #a,b,...,z,aa,ab,...,az
    for i in range(702, 728):
        print labelP(i, alphabet)
        #aaa,aab,...,aaz
    for i in range(18278, 18304):
        print labelP(i, alphabet)
        #aaaa, aaab, ... , aaaz

