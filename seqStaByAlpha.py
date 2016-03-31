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

def readGenome(genome):
    aDict = {}
    key = ''
    for line in open(genome):
        if line[0] == '>':
            if key:
                aDict[key] = ''.join(alist)
            #--------------------------------
            key = line[1:-1]
            alist = []
        else:
            alist.append(line.upper().strip())
    #-------------------------------------------
    aDict[key] = ''.join(alist)
    return aDict
#--------------------------------

def getSeq(bed,  aDict):
    '''
    bed: 0-based,  half-open
    '''
    output = bed + '.fa'
    fh = open(output, 'w')
    for line in open(bed):
        lineL = line.strip().split('\t')
        chr, start, end, name = lineL[:4] 
        start = int(start)
        end   = int(end)
        print >>fh, '>%s\n%s' % (name,  aDict[chr][start:end])
    fh.close()
#-------------------------------------

def countCpG(bed):
    '''
    CG count
    >5 CpG per 100bp
    1-5 CpG per 100bp
    <1 CpG per 100bp
    '''
    aDict = {'Regions with < 1 CpG per 100bp': 0, 
             'Regions with 1-5 CpG per 100bp': 0, 
             'Regions with > 5 CpG per 100bp': 0}
    for line in open(bed+'.fa'):
        if line[0] != '>':
            line = line.strip()
            lenline = len(line)
            CpG = line.count('CG')
            count = CpG * 100.0 / lenline
            if CpG * 100.0 / lenline < 1:
                aDict['Regions with < 1 CpG per 100bp'] += 1
            elif count <= 5:
                aDict['Regions with 1-5 CpG per 100bp'] += 1
            else:
                aDict['Regions with > 5 CpG per 100bp'] += 1 
        #===============================================
    #--------------------------------------------------------
    fh = open(bed + '.cpg',  'w')
    print >>fh,  "Sample\tRegions with < 1 CpG per 100bp\tRegions with\
 1-5 CpG per 100bp\tRegions with > 5 CpG per 100bp"
    print >>fh, "%s\t%d\t%d\t%d" % (bed, \
            aDict['Regions with < 1 CpG per 100bp'], \
            aDict['Regions with 1-5 CpG per 100bp'], \
            aDict['Regions with > 5 CpG per 100bp'])
    fh.close()
#------------------------------------------
def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) < 3:
        print >>sys.stderr, 'Using python %s genome \
pos_bed1[,pos_bed2...]' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------------------
    genome = sys.argv[1]
    bedL = sys.argv[2:]
    aDict = readGenome(genome)
    for bed in bedL:
        getSeq(bed,  aDict)
        countCpG(bed)
if __name__ == '__main__':
    main()

