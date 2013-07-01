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

This is used to compute the co-existence of eithe two motifs in a
file.

Input file format:

	FASTA ID        Offset  Sequence        Motif Name      Strand  MotifScore
	P_10753_1_2_3   -2      GGACACTG        1-GGACTCDV      +   6.214134
	P_22195 -39     GGACTGTG        1-GGACTCDV      +       6.513756
	P_22195 2       GAACTGCA        1-GGACTCDV      +       4.727761
	P_22195 10      GGACAGGC        1-GGACTCDV      +       5.319863
	P_22195 64      GGACTCGC        1-GGACTCDV      +       7.140175
	P_33660_1_2_3_4_5_6     -2      GGACTCGG        1-GGACTCDV      +   	7.163865
	P_19079 -4      GGACATCA        1-GGACTCDV      +       4.602924
	P_29105 -59     GGACTGGG        1-GGACTCDV      +       6.403520
	P_29105 -26     GGACAGGT        1-GGACTCDV      +       4.533459
	P_32005_1_2     -2      GAACTGGA        1-GGACTCDV      +   	5.355430
	P_32005_1_2     61      GGACAGTA        1-GGACTCDV      +   	5.377363
	P_21555 -2      GGACTGTG        1-GGACTCDV      +       6.513756
	P_21555 39      GGACTTGG        1-GGACTCDV      +       6.366986
	P_25749_1_2     -2      GGACTCGG        1-GGACTCDV      +   	7.163865
	P_25749_1_2     16      GGACAGAG        1-GGACTCDV      +   	5.443637

Output file format:
    Per motif1  motif2  motif3  ....
motif1  100%    90% 80% ...
motif2  80% 100%    90% ...
motif3  70% 60% 100%    ...
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
        metavar="FILEIN", help="The output of findMotifs.pl or \
findMotifsGenome.pl.")
    parser.add_option("-m", "--max-distance", dest="dist",
        metavar="FILEIN", default=100000, help="The maximum distance \
to define co-existence. If the distance between two motifs are larger \
than the one given here, they will not be taken as co-linear. Default \
100000 (a meaningless large value). ")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-n", "--column_number_for_name", dest="name_col",
        default=1, help="The column number which represent \
the name of each region. Default 1. (1-based)")
    parser.add_option("-p", "--column_number_for_pos", dest="pos_col",
        default=2, help="The column number which represent \
the position of each motif. Default 2. (1-based)")
    parser.add_option("-q", "--column_number_for_motif", dest="motif_col",
        default=4, help="The column number which represent \
the name of each motif. Default 4. (1-based)")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    dist = int(options.dist)
    verbose = options.verbose
    debug = options.debug
    name_col = int(options.name_col) -1
    pos_col = int(options.pos_col) -1 
    motif_col = int(options.motif_col) -1
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    motifCntD = {}
    #aDict = {region_name:{motif1:[pos1, pos2], \
    #    motif2:[pos]}}
    header = 1
    for line in fh:
        if header:
            header -= 1
            continue
        #-----------------------------
        lineL = line.split()
        name = lineL[name_col]
        pos  = lineL[pos_col]
        motif = lineL[motif_col]
        if name not in aDict:
            aDict[name] = {}
        if motif not in aDict[name]:
            aDict[name][motif] = []
        aDict[name][motif].append(int(pos))
        if motif not in motifCntD:
            motifCntD[motif] = set()
        motifCntD[motif].add(name)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    motifKeyL = motifCntD.keys()
    motifKeyL.sort()
    len_motifkeyL = len(motifKeyL)
    combinedKeyL = []
    for i in range(len_motifkeyL-1):
        for j in range(i+1, len_motifkeyL):
            combinedKeyL.append('^'.join([motifKeyL[i], motifKeyL[j]]))
    for i in motifKeyL:
        combinedKeyL.append('^'.join([i, i]))
    #---------------------------------------------------
    combinedKeyD = {}
    for i in combinedKeyL:
        combinedKeyD[i] = 0
    for key, valueS in motifCntD.items():
        motifCntD[key] = len(valueS)
    #print >>sys.stderr, motifCntD
    for name, motifD in aDict.items():
        for ck in combinedKeyL:
            twomotif = 2
            ck1, ck2 = ck.split('^')
            if ck1 not in motifD or ck2 not in motifD:
                continue
            ck1PosL = motifD[ck1]
            lenck1 = len(ck1PosL)
            ck2PosL = motifD[ck2]
            if ck1 == ck2:
                if lenck1 > 1:
                    for i in range(lenck1-1):
                        if abs(ck1PosL[i] - ck1PosL[i+1]) <= dist:
                            combinedKeyD[ck] += 1
                            break
                #-----------------------------------------
            else:
                yes = 0
                for pre in ck1PosL:
                    for post in ck2PosL:
                        if abs(pre-post) <= dist:
                            combinedKeyD[ck] += 1
                            yes = 1
                            break
                    #--------------------------
                    if yes:
                        break
                #--------------------------
            #--end one motif combination----------------
        #------end one region-------------------------------
    #--end all regions----------------------------
    if verbose:
        print >>sys.stderr, combinedKeyD
        print >>sys.stderr, motifCntD
    print 'Per\t%s' % '\t'.join(motifKeyL)
    for i in range(len_motifkeyL):
        ck1 = motifKeyL[i]
        tmpL = [ck1]
        for j in range(len_motifkeyL):
            ck2 = motifKeyL[j]
            if i<= j:
                ck = ck1+'^'+ck2
            else:
                ck = ck2+'^'+ck1
            tmpL.append(str(combinedKeyD[ck]*1.0/motifCntD[ck1]))
        #------------------------------------
        print '\t'.join(tmpL)
    #-------------------------------------
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



