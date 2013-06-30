#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2013, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
'''
Functionla description

Input fie format (Duplicate names are allowed, the one with the
largest value will be kept.)
NM_001177713_7331.UTR5.5        0
NM_001177713_7331.UTR5.4        0
NM_001177713_7331.UTR5.3        0
NM_001177713_7331.UTR5.2        0
NM_001177713_7331.UTR5.1        0
NM_001177713_7331.Coding_exon.50        0
NM_001177713_7331.Coding_exon.49        0
NM_001177713_7331.Coding_exon.48        0
NM_001177713_7331.Coding_exon.47        0
.
.
.
.
NM_001177713_7331.Coding_exon.4 0
NM_001177713_7331.Coding_exon.3 0
NM_001177713_7331.Coding_exon.2 0
NM_001177713_7331.Coding_exon.1 0
NM_001177713_7331.UTR3.20       0
NM_001177713_7331.UTR3.19       0
NM_001177713_7331.UTR3.18       0
NM_001177713_7331.UTR3.17       0
NM_001177713_7331.UTR3.16       0
.
.
.
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    desc = ""
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="")
    parser.add_option("-x", "--x-label", dest="x_label",
        metavar="1", default='1', help="The part of strings used for x-label. \
Starts with 1")
    parser.add_option("-y", "--y-label", dest="y_label",
        metavar="2,3", default='2,3', help="The part of strings used for y-label. \
Starts with 1")
    parser.add_option("-s", "--separtor", dest="separtor",
        metavar=".", default='.',  help="String separtor.")
    parser.add_option("-k", "--key_y_label", dest="key",
        metavar="UTR5.20.Coding_exon.50.UTR3.20",
        default='UTR5.20.Coding_exon.50.UTR3.20',  
        help="UTR5, Coding_exon, UTR3 represents label type. It can be \
any given string. 20,50,20 means the numer of each type before them. \
In this example , we will get <UTR5.1,UTR5.2,\
...UTR5.20,Coding_exon.1,...,Coding_exon.50,UTR3.1,UTR3.2...,UTR3.20>\
for y_label.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
#aDict = {gene:{'UTR5.1':1, }}

def output(aDict, keyL):
    print 'gene\t%s' % '\t'.join(keyL)
    geneL = aDict.keys()
    geneL.sort()
    for gene in geneL:
        tmpL = [aDict[gene].get(i,'0') for i in keyL]
        print '%s\t%s' % (gene, '\t'.join(tmpL))
#------------------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    verbose = options.verbose
    debug = options.debug
    sep = options.separtor
    key = options.key.split('.') #UTR5.20.Coding_exon.50.UTR3.20
    lenkey = len(key)
    keyL = []
    for i in range(1,lenkey,2):
        tmpL = [sep.join([key[i-1],str(j+1)]) for j in range(0,int(key[i]))]
        keyL.extend(tmpL)
    #print keyL
    #sys.exit()
    x_label = [int(i)-1 for i in options.x_label.split(',')]
    y_label = [int(i)-1 for i in options.y_label.split(',')]
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    #oldgene = ''
    #keyL = []
    
    for line in fh:
        name, value = line.split()
        nameL = name.split(sep)
        gene = sep.join([nameL[i] for i in x_label])
        type = sep.join([nameL[i] for i in y_label])
        assert type in keyL, type

        #if oldgene == "":
        #    oldgene = gene
        #if oldgene == gene:
        #    if type not in keyL:
        #        keyL.append(type)
        #    else:
        #        assert type == keyL[-1]
        if gene not in aDict:
            aDict[gene] = {}
        #if type not in aDict[gene]:
        if type not in aDict[gene]:
            aDict[gene][type] = value
        elif value > aDict[gene][type]:
            aDict[gene][type] = value
        #else:
        #    print >>sys.stderr, "Duplicate type name", line
        #    sys.exit()
        #--------------------------------------------
    #-------------END reading file----------
    #print >>sys.stderr, aDict['NM_001177713_7331']
    output(aDict, keyL)
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    if verbose:
        print >>sys.stderr,\
            "--Successful %s" % strftime(timeformat, localtime())
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()



