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
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

desc='''
This is used to get fasta sequence based on given locus.
'''

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A FASTA file. Sequences in multiple \
lines are supported.")
    parser.add_option("-l", "--locus", dest="locus",
        metavar="LOCUS", help="A file with given locus one each line.")
    parser.add_option("-s", "--sep", dest="sep",
        metavar="SEP",default=" ", help="The separtor for the name \
FASTA seq in given FASTA file. Default one space.")
    parser.add_option("-o", "--op", dest="op",
        metavar="OP",default="include", help="Default <include> which \
means only those sequences with a name matching to locus will be \
outputted. <exclude> means sequence with locus not in locus-file will \
be output.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readFasta(file, sep, locusD, op):
    '''
    '''
    #------------------------------------------
    aDict = {}
    for line in open(file):
        if line[0] == '>':
            saveThis = 1
            locus = (line[1:].split(sep)[0]).strip()
            if op=="include" and locus not in locusD:
                saveThis = 0
            elif op=='exclude' and locus in locusD:
                saveThis = 0
            else:
                aDict[locus] = ''
        elif saveThis:
            aDict[locus] += line.strip() 
        #-------------------------------------------
    #-------------------------------------------
    return aDict
#-----------------------------------------------------------


def main():
    #--------------------------------------
    options, args = cmdparameter(sys.argv)
    file = options.filein
    locus = options.locus
    sep = options.sep
    op = options.op
    locusD = dict([(line.strip(), 1) for line in open(locus)])
    seqDict = readFasta(file, sep, locusD, op)
    keyL = seqDict.keys()
    keyL.sort()
    for key in keyL:
        print '>%s\n%s' % (key, seqDict[key])

if __name__ == '__main__':
    main()




