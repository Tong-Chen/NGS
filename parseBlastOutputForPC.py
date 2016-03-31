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

def readBlast(file):
    aDict = {}
    for line in open(file):
        lineL = line.split('\t')
        gi, identity = lineL[1], float(lineL[2])
        gi = gi.split('|', 2)[1]
        assert gi not in aDict
        aDict[gi] = identity
    return aDict
#------------------------------------
def readDbcmd(file):
    aDict = {}
    locus = ''
    for line in open(file):
        if line[0] == '>':
            if locus:
                aDict[locus] = ''.join(aDict[locus])
            locus = line.split('|', 2)[1]
            aDict[locus] = []
        else:
            aDict[locus].append(line.strip())
    #------------------------------------
    #----last line--------------
    if locus:
        aDict[locus] = ''.join(aDict[locus])
    return aDict
#---------------------------------------
def readSql(file,  head=1):
    aDict = {}
    for line in open(file):
        if head:
            head -= 1
            continue
        #----------------------------
        value, key = line.split()
        if key not in aDict:
            aDict[key] = [value]
        else:
            aDict[key].append(value)
    #-------------------------------------------
    return aDict
#-----------------------------------------------
def main():
    if len(sys.argv) != 4:
        print >>sys.stderr, "Preprocess blast hit for following analysis."
        print >>sys.stderr, 'Using python %s blastoutput \
blastdbcmdoutput[gi-fasta-seq] mysqloutput[gi-taxid]' % sys.argv[0]
        sys.exit(0)
    #---------------------------------------------------
    blastDict = readBlast(sys.argv[1])
    dbcmdDict = readDbcmd(sys.argv[2])
    sqlDict   = readSql(sys.argv[3])
    #-----------------------------------
    for key, valueL in sqlDict.items():
        if len(valueL) > 1:
            tmpgi = valueL[0]
            tmpvalue = blastDict[tmpgi]
            for gi in valueL:
                curvalue = blastDict[gi]
                if curvalue > tmpvalue:
                    tmpgi = gi
                    tmpvalue = curvalue 
            #------------------------
            sqlDict[key] = tmpgi
        else:
            sqlDict[key] = sqlDict[key][0]
    #----------------------------------
    for key, value in sqlDict.items():
        if value in dbcmdDict:
            print ">%s\n%s" % (key, dbcmdDict[value])
        else:
            print "Wrong %s" % value
            sys.exit(1)
#----------------------------------------------
if __name__ == '__main__':
    main()

