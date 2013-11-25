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
This is used to get hairpin sequence for given species.
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
        metavar="FILEIN", help="A FASTA file.")
    parser.add_option("-s", "--species", dest="species",
        metavar="species", help="The short name of a species. \
Like mmu.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readFasta(file, spe):
    '''
    '''
    #------------------------------------------
    aDict = {}
    saveThis = 0
    start = '>' + spe + '-'
    for line in open(file):
        if line[0] == '>' and not line.startswith(start):
            saveThis = 0
        elif line[0] == '>' and line.startswith(start):
            saveThis = 1
            locus = line[1:].strip()
            aDict[locus] = ''
        elif saveThis:
            aDict[locus] += line.strip() 
        #-------------------------------------------
    #-------------------------------------------
    for key,item in aDict.items():
        print '>%s\n%s' % (key, item)
#-----------------------------------------------------------


def main():
    #--------------------------------------
    options, args = cmdparameter(sys.argv)
    file = options.filein
    spe = options.species
    readFasta(file, spe)

if __name__ == '__main__':
    main()




