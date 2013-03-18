#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division#, with_statement
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
import sys


def countCpG(file):
    '''
    CG count
    >5 CpG per 100bp
    1-5 CpG per 100bp
    <1 CpG per 100bp
    '''
    aDict = {'Regions with < 1 CpG per 100bp': 0, 
             'Regions with 1-5 CpG per 100bp': 0, 
             'Regions with > 5 CpG per 100bp': 0}
    bDict = {'Regions with < 1 CpG per 100bp': [0,0], 
             'Regions with 1-5 CpG per 100bp': [0,0], 
             'Regions with > 5 CpG per 100bp': [0,0]}
    for line in open(file):
        if line[0] != '>':
            line = line.strip()
            lenline = len(line)
#            CpG1 = line.count('CG')
#            CpG2 = line.replace('CG','  ').count('GC')
            c = line.count('C') * 1.0 / lenline
#            print c
            g = line.count('G') * 1.0 / lenline
            CpG = line.count('CG') #+ line.replace('CG','  ').count('GC')
            count = CpG * 100.0 / lenline
            if CpG * 100.0 / lenline < 1:
                aDict['Regions with < 1 CpG per 100bp'] += 1
                bDict['Regions with < 1 CpG per 100bp'][0] += c
                bDict['Regions with < 1 CpG per 100bp'][1] += g
            elif count <= 5:
                aDict['Regions with 1-5 CpG per 100bp'] += 1
                bDict['Regions with 1-5 CpG per 100bp'][0] += c
                bDict['Regions with 1-5 CpG per 100bp'][1] += g
            else:
                aDict['Regions with > 5 CpG per 100bp'] += 1 
                bDict['Regions with > 5 CpG per 100bp'][0] += c
                bDict['Regions with > 5 CpG per 100bp'][1] += g
        #===============================================
    #--------------------------------------------------------
    fh = open(file + '.cpg',  'w')
    print >>fh,  "Sample\tRegions with < 1 CpG per 100bp\tRegions with\
 1-5 CpG per 100bp\tRegions with > 5 CpG per 100bp"
    print >>fh, "%s\t%d\t%d\t%d" % (file, \
            aDict['Regions with < 1 CpG per 100bp'], \
            aDict['Regions with 1-5 CpG per 100bp'], \
            aDict['Regions with > 5 CpG per 100bp'])
    #----------------------------------------------------------
    print >>fh, "%s\t%f\t%f\t%f" % (file, \
            bDict['Regions with < 1 CpG per 100bp'][0]/aDict['Regions with < 1 CpG per 100bp'], \
            bDict['Regions with 1-5 CpG per 100bp'][0]/aDict['Regions with 1-5 CpG per 100bp'], \
            bDict['Regions with > 5 CpG per 100bp'][0]/aDict['Regions with > 5 CpG per 100bp'])
    #----------------------------------------------------------
    print >>fh, "%s\t%f\t%f\t%f" % (file, \
            bDict['Regions with < 1 CpG per 100bp'][1]/aDict['Regions with < 1 CpG per 100bp'], \
            bDict['Regions with 1-5 CpG per 100bp'][1]/aDict['Regions with 1-5 CpG per 100bp'], \
            bDict['Regions with > 5 CpG per 100bp'][1]/aDict['Regions with > 5 CpG per 100bp'])

    fh.close()
#------------------------------------------
def main():
    print >>sys.stderr, "Count the number of CG in the specified regions."
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) < 3:
        print >>sys.stderr, 'Using python %s  \
fasta_file1[, fasta_file2...]' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------------------
    bedL = sys.argv[1:]
    for bed in bedL:
        countCpG(bed)
if __name__ == '__main__':
    main()

