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
    parser.add_option("-n", "--number-bins", dest="nBins",
        metavar="20,50,20", default="20,50,20", 
        help="The numbers for UTR5, CDS, UTR3 \
regions")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def getBinsType(innerD, type, nBin, gene):
    typeL = innerD[type]
    type_len = sum([int(i[2])-int(i[1]) for i in typeL])
    size = type_len / nBin
    if size == 0:
        size = 1
    cnt = 0
    for i in typeL:
        #print i[1], i[2], size
        for j in range(int(i[1]), int(i[2]), size):      
            cnt += 1
            end = j + size
            #print j, end
            if end + size <= int(i[2]):
                tmpL = [i[0], str(j), str(end), \
                    gene+'.'+type+'.'+str(cnt), i[4], i[5]]
                print '\t'.join(tmpL)
            else:
                tmpL = [i[0], str(j), i[2], \
                    gene+'.'+type+'.'+str(cnt), i[4], i[5]]
                print '\t'.join(tmpL)
                #cnt += 1
                break
                #------------------------------
            #print '\t'.join(tmpL)
            #cnt += 1
    #--------------------------------------------------
    for i in range (cnt, nBin):
        lastCnt = int(tmpL[3].split('.')[2])
        lastCnt += 1
        tmpL[3] = gene+'.'+type+'.'+str(lastCnt)
        print '\t'.join(tmpL)


#---------------------------------------

def getBins(aDict, nBins):
    for gene, innerD in aDict.items():
        if 'UTR5' in innerD:
            getBinsType(innerD, 'UTR5', nBins[0], gene)
        if 'Coding_exon' in innerD:
            getBinsType(innerD, 'Coding_exon', nBins[1], gene)
        if 'UTR3' in innerD:
            getBinsType(innerD, 'UTR3', nBins[2], gene)
#-----------------------------------------------------------
def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    nBins = [int(i) for i in options.nBins.split(',')]
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    for line in fh:
        lineL = line.split()
        keyL = lineL[3].split('.')
        gene = keyL[0]
        type = keyL[1]
        if gene not in aDict:
            aDict[gene] = {}
        #---------------------------------
        if type in aDict[gene]:
            aDict[gene][type].append(lineL)
        else:
            aDict[gene][type] = [lineL]
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    #print aDict['NM_027855_3']
    getBins(aDict, nBins)
    #---------------------------------------
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



