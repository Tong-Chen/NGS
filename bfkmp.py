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

def index(subseq, seq):  
    ''''' 
    index(subseq, seq) -->a list of numbers or -1 
    Return an index of [subseq] in the [seq]. 
    Or -1 if [subseq] is not a subsequence of [seq]. 
    The time complexity of the algorithm is O(n*m), where 
    n, m = len(seq), len(subseq). 
    >>>index('12', '0112') 
    [2] 
    >>>index([1,2], [011212]) 
    [2, 4] 
    >>>index('13', '0112') 
    -1 
    '''  
    i, n, m = -1, len(seq), len(subseq)  
    index = []  
    try:  
        while True:  
            i = seq.index(subseq[0], i+1, n - m + 1)  
            if subseq == seq[i:i+m]:  
                index.append(i)  
    except ValueError:  
        return index if len(index) > 0 else -1  
def subseqInSeq(subseq, seq):  
    ''''' 
    subseqInSeq(subseq, seq) ---> list or -1 
    The same as index. 
    '''  
    indexList = []  
    m = len(subseq)  
    subseqRepla = '*' * m  
    while subseq[0] in seq:  
        index = seq.index(subseq[0])  
        if subseq == seq[index:index+m]:  
            indexList.append(index)  
            seq = seq.replace(subseq, subseqRepla, 1)  
        else:  
            seq = seq.replace(subseq[0], '*', 1)  
    return (indexList if len(indexList) > 0 else -1)   
def main():  
    print index('ab', 'abcdab')  
    print subseqInSeq('ab', 'abcdab')  
