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
        metavar=".", default='.',  help="The part of strings used for x-label. \
Starts with 1")
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
    x_label = [int(i)-1 for i in options.x_label.split(',')]
    y_label = [int(i)-1 for i in options.y_label.split(',')]
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    oldgene = ''
    keyL = []
    for line in fh:
        name, value = line.split()
        nameL = name.split(sep)
        gene = sep.join([nameL[i] for i in x_label])
        type = sep.join([nameL[i] for i in y_label])
        if oldgene == "":
            oldgene = gene
        if oldgene == gene:
            if type not in keyL:
                keyL.append(type)
            else:
                assert type == keyL[-1]
        if gene not in aDict:
            aDict[gene] = {}
        #if type not in aDict[gene]:
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



