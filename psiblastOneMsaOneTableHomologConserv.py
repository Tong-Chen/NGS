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

def readsimple(hitfile):
    '''
    hitfile means the files producted by psiblastOneMsaOneTable.py.
    aDict = {searchlocus:set(hitlocus)}
    '''
    aDict = {}
    for line in open(hitfile):
        if line[0] == '=':
            slocus = ((line[1:].split()[0]).rsplit('.', 1))[0]
            num = 2 #To skip the following two lines which are
                    #search seq
            if slocus not in aDict:
                aDict[slocus] = set()
        elif num: #skip two lines
            num -= 1
        else: #only save the locus
            ###Here is the changing part--
            ###To decide which part of locus you want
            if line[0] == '>':
                aDict[slocus].add(line[1:-1])
        #-------------------------------------------
    #--------------END for------------------           
    return aDict
#---------------------------END file-------------------

def output(hitDict, deleteself = 0):
    '''
    hitDict = {searchlocus: set(hitlocus)}
    '''
    for key, valueS in hitDict.items():
        lenvalueS = len(valueS)
        if lenvalueS == 0:
            print '>%s* %d' % (key, lenvalueS)
        else:
            if key not in valueS:
                print '>%s* %d' % (key, lenvalueS)
                print '#'.join(valueS)
            else:
                valueS.remove(key)
                newlenvalues = len(valueS)
                assert (newlenvalues == lenvalueS-1), 'Wrong operation'
                print '>%s %d' % (key, newlenvalues)
                if newlenvalues:
                    print '#'.join(valueS)
        #--------------------------------------
    #-------------------------------------------------               
#---------------------------END output-------------------
def testADict(aDict):
    for key, value in aDict.items():
        print '>%s' % key
        print '%s' % '#'.join(value)

def main():
    print >>sys.stderr, "Print the result to screen"
    print >>sys.stderr, "This used to parse the result produced by \
psiblastOneMsaOneTable.py. Specially this one is for ortholog \
conservation analysis. "
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #--------------------------------------
    aDict = readsimple(sys.argv[1])
    output(aDict)


if __name__ == '__main__':
    main()

