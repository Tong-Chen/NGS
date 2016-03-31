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

def gi2taxonID(file):
    #file = '/home/CT/server/ncbi/gi_taxid_prot.dmp'
    aDict = {}
    for line in open(file):
        gi, taxonID = line.strip().split()
        aDict[gi] = taxonID
    #------------------------
    return aDict
#-----------------------------------------------------

def taxonID2Name(file= '/home/CT/server/ncbi/names.dmp'):
    aDict = {}
    ##only take scientific name----
    for line in open(file):
        if line.find('scientific name') != -1:
            taxonID, taxonName, nouse = line.split('\t|\t',2)
            aDict[taxonID] = taxonName
    #------------------------
    return aDict
#--------------------------------------------------------

def main():
    #id2name = taxonID2Name()
    #for key, value in id2name.items():
    #    print "%s\t%s" % (key, value)
    gi2id = gi2taxonID()
    for key, value in gi2id.items():
        print "%s\t%s" % (key, value)
if __name__ == '__main__':
    main()

