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

Input files:

    chr1	23853198	23853221	1	CTGCAAAGGCATCTAGTTCCCAT	-
    chr1	34872326	34872349	2	AGCACCAGCAAGCACACCCTGCA	-
    chr1	37932703	37932720	3	TCGGCAGCAAATTCTGCCCCG	+
    chr1	40829459	40829463	3	TCGGCAGCAAATTCTGCCCCG	+



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
        metavar="FILEIN", help="A standard bed file with \
six columns. Multiple separated regions with same names will \
be preocessed together as you would expect. \
Strand information will be considered. \
Extra columns will be ignored.")
    parser.add_option("-n", "--number-bins", dest="nBins",
        metavar="NUMBER-BINS", default="20", 
        help="Set the number of bins for given regions")
    parser.add_option("-s", "--strand", dest="strand",
        default=1, help="Consider strand information or not. \
Default TRUE. Strings or numbers that represent FALSE can turn \
off this parameter.")
    parser.add_option("-u", "--upstream", dest="up",
        default=2000, help="The number of nucleotides extended \
to upstream. Default 2000 for CpG shore.")
    parser.add_option("-d", "--downstream", dest="dw",
        default=2000, help="The number of nucleotides extended \
to downstream. Default 2000 for CpG shore.")
    parser.add_option("-N", "--number-bins-flank", dest="nBins_flank",
        metavar="NUMBER-BINS-FLANK", default="40", 
        help="Set the number of bins for flanking regions")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def getBinsType(typeL, nBin, gene, strand):
    #add sort operation
    #comment out sort here
    #typeL.sort(key=lambda x:int(x[1]))
    #print >>sys.stderr, typeL
    type_len = sum([int(i[2])-int(i[1]) for i in typeL])
    size = type_len / nBin
    mod = type_len % nBin
    if mod:
        size = size + 1
    #---------------------------
    #if size < 2:
    #if size == 1:
    #    size = 2
    #elif size == 0:
    #    size = 1
    if strand == '+':
        cnt = 0
        time = 1
    elif strand == '-':
        cnt = nBin + 1
        time = -1
    else:
        print >>sys.stderr, "Wrong strand type %s" % strand
    remain = 0
    lenTypeL = len(typeL)
    for index in range(lenTypeL):
        i = typeL[index]
        start = int(i[1])
        while True:
            if remain == 0:
                cnt += 1 * time
                end = start + size
            else:
                end = start + remain
                remain = 0
            #if index != lenTypeL -1: 
            if end  <= int(i[2]):
                tmpL = [i[0], str(start), str(end), \
                    gene+'.'+str(cnt)]
                tmpL.extend(i[4:])
                print '\t'.join(tmpL)
                start = end
            else:
                remain = end - int(i[2])
                tmpL = [i[0], str(start), i[2], \
                    gene+'.'+str(cnt)]
                tmpL.extend(i[4:])
                print '\t'.join(tmpL)
                start = int(i[2])
                break
            #------------------------------------------------
            if start >= int(i[2]):
                break
        #--------------------------------------------------
    #--------------------------------------------------
    if strand == '+':
        for i in range (cnt, nBin):
            lastCnt = int(tmpL[3].split('.')[-1])
            lastCnt += 1
            tmpL[3] = gene+'.'+str(lastCnt)
            print '\t'.join(tmpL)
    else:
        for i in range (cnt, 1, -1):
            lastCnt = int(tmpL[3].split('.')[-1])
            lastCnt -= 1
            assert lastCnt >= 1, "%s\t%d\t%d" % (type, cnt, lastCnt)
            tmpL[3] = gene+'.'+str(lastCnt)
            print '\t'.join(tmpL)
    #-----------------------------------------------------

#---------------------------------------

def getBins(aDict, nBins, strandD, nameAdd=''):
    for gene, innerL in aDict.items():
        #if not strandD:
        #    getBinsType(innerL, nBins, gene+nameAdd, '+')
        #else:
        getBinsType(innerL, nBins, gene+nameAdd, strandD.get(gene, '+'))
#-----------------------------------------------------------
def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    nBins = int(options.nBins)
    up = int(options.up)
    dw = int(options.dw)
    nBins_flank = int(options.nBins_flank)
    verbose = options.verbose
    debug = options.debug
    strand = int(options.strand)
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    strandD = {}
    for line in fh:
        lineL = line.split()
        gene = lineL[3]
        #keyL = lineL[3].split('.')
        #gene = keyL[0]
        #type = keyL[1]
        if gene not in aDict:
            aDict[gene] = [lineL]
        else:
            aDict[gene].append(lineL)
        #---------------------------------
        if strand and len(lineL) >= 6:
            strandD[gene] = lineL[5]
        #if type in aDict[gene]:
        #    aDict[gene][type].append(lineL)
        #else:
        #    aDict[gene][type] = [lineL]
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    # sort aDict
    for gene, typeL in aDict.items():
        typeL.sort(key=lambda x:int(x[1]))
        #type_len = sum([int(i[2])-int(i[1]) for i in typeL])
    #print aDict['NM_027855_3']
    #print aDict
    if up or dw:
        getBins(aDict, nBins, strandD, "@given")
    else:
        getBins(aDict, nBins, strandD)
    #---------------------------------
    if up:
        upD = {}
        for gene, typeL in aDict.items():
            #upD[gene] = [typeL[0][0]]
            newend = int(typeL[0][1])
            newstart = newend - up
            if newstart < 0:
                newstart = 0
            #upD[gene].end(newstart, newend)
            upD[gene] = [[typeL[0][0], str(newstart), str(newend)]]
            upD[gene][0].extend(typeL[0][3:])
        getBins(upD, nBins_flank, strandD, "@up")
    #---------------------------------------
    if dw:
        dwD = {}
        for gene, typeL in aDict.items():
            #dwD[gene] = [typeL[-1][0]]
            newstart = int(typeL[-1][2])
            newend = newstart + dw
            #if newstart < 0:
            #    newstart = 0
            #dwD[gene].append(newstart, newend)
            #dwD[gene].extend(typeL[-1][3:])
            dwD[gene] = [[typeL[-1][0], str(newstart), str(newend)]]
            dwD[gene][0].extend(typeL[-1][3:])
        getBins(dwD, nBins_flank, strandD, "@dw")
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



