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
This is written to transfer BAM file to wig file. 
Currently this script only focus on Pair-end RNA-Seq reads.
One feature of this program is that it can recover the coverage for
<insert regions> between two paired reads.

For example,  we have a properly mapped reads showed below


==============================================  Genome
^^^^^^                                  $$$$$$  Two reads 
      ++++++++++++++++++++++++++++++++++        Insert regions

Traditionally, the output wig file for regions labeld '+' would be 0.
Actually,  these regions are also covered by reads. Thus,  using this
program,  you can get <true> coverage for <insert regions>.


==============================================  Genome
^^^^^^         ^^^^^^                           Two reads 
      +++++          +++++                      Extend regions


It was also planned to transfer BAM with single-end reads to wig.
However, this is not finished yet. But other substitutaion tools or 
combined tools are available in this directory can deal with this type
of transferation.  

1.PE extend-unextend strand-unstrand finished
2.SE extend-unextend strand-unstrand finished
  a. 2013-12-26 SE unstranded RNA extend make through the test
STH unexpected:
1.Mapped reads with flag 83 even when two reads mapped to different
chromosomes.

'''

import collections
import sys
import re
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
DEBUG = 0

def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    desc = "This is written to transfer BAM file to \
wig file. It can deal with single-end bam and pair-end bam and \
strand-specific bam and the combination of them. For pair-end bam, a \
GTF file is needed for RNA-Seq. The program will try to assign reads for middle \
regions."
    usages = "%prog [-i SAM file] -g GTF -n RNA -t PE -l fr-unstranded -e 1"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A SAM file with or without header,  or \
- means STDIN. No header lines allowed and file must be sorted by chromosome")
    parser.add_option("-g", "--gtf", dest="gtf",
        metavar="GTF", help="When -t is PE and -e is positive, this \
GTF should be supplied. It canbe standard GTF downlaoded from UCSC. \
However, the one outputted by you RNA-Seq would be better (Only with \
expressed transcripts is preferred).")
    parser.add_option("-n", "--nucleotide-type", dest="nt",
        metavar="RNA/DNA", default="RNA", 
        help="Default RNA. DNA means ChIP-Seq, RNA means RNA-Seq. \
Currently only RNA is supported.")
    parser.add_option("-t", "--seq-type", dest="seq_Type",
        metavar="PE/SE", default="PE", 
        help="Default PE. PE for pair-end reads and SE for single-end reads.")
    parser.add_option("-l", "--library-type", dest="lt",
        default="fr-unstranded", 
        help="Default fr-unstranded. fr-unstranded,fr-firststrand,\
fr-secondstrand. In default value, all exons in GTF will be treated as \
positive strand exons.")
    parser.add_option("-e", "--extend", dest="extend", default=0,
        type='int', help="Default 0 means no extend. A positive number \
indicates the expected legth should be given for SE data. \
For PE reads, any positive number would be enough to let \
the progeam to fill in the blank between two PE reads assisted by \
GTF. 0 means no extending or filling. ")
    parser.add_option("-c", "--chrom-size", dest="chromSize",
        help="Only needed when <nucleotide-type> is DNA. \
If -t is SE and -e larger than 0 and no header in SAM \
file. This should be given. One can use \
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
'select chrom, size from mm10.chromInfo' > mm10.genome \
to extract chromosome size.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information.")
    parser.add_option("-d", "--debug", dest="debug",
        default=0, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def computeRegion(start, cigarL, name): 
    # cigarL:  [('20','M'),('1000','N'),('80','M')]
    regionL = []
    append = regionL.append
    len = 0
    for i in cigarL:
        if i[1] == 'N': #begin settle counts
            assert len > 0, "Unexpected cigarL %s " % name
            append([start, start+len])
            start = start + len + int(i[0])
            len = 0
        else:
            if i[1] == 'M' or i[1] == 'D':
                len += int(i[0])
            elif i[1] == 'I':
                pass # no adding length
            else:
                print >>sys.stdout, "Unconsidered cigars %s" % name
        #------------------------------------
    #--------------END for-------------------
    #------If no 'N' or deal with the part after the last 'N'--
    append([start, start+len])
    return regionL
    #--------------------------------------------------------------------
#--------END computeRegion-------------------------------------------------

def computeWigDict(wigDict, pairL, extend):
    pairL = pairL
    #wigDict = wigDict
    #wigDict = {pos:{+:[+,+_e], '-':[ -,-_e]}}
    #pairL   = [[chr,flag,[[start,end],...],xs], [chr,flag,[[start,end],...],xs]]
    #---Get relative position for two reads--------------
    # -1 means reads with smaller coordinate (more left  in Genome)
    #  0 means reads with bigger  coordinate (more right in Genome)
    # These two numbers are also used for sorting two reads and
    # extract neighbor mapped regions of two reads.
    for i in pairL:
        i[1] = -1 if i[1] & 0x20 else 0
        #if i[1] & 0x20: #mate reverse 99(pPR1),163(pPR2)
        #    i[1] = -1 
        #elif i[1] & 0x10: #self reverse 147(pPr2), 83(pPr1)
        #                 #reads that self reversed located right in
        #                 #Genome
        #    i[1] = 0
    #----------------------------------------------------
    #----------------------------------------------------------- 
    assert pairL[0][0] == pairL[1][0], pairL
    assert pairL[0][3] == pairL[1][3], pairL
    xs = pairL[0][3]
    pairL.sort(key=lambda i:i[1])
    #print >>sys.stderr, pairL
    leftReadsMaxCor = pairL[0][2][-1]
    rightReadsMinCor = pairL[1][2][0]
    #--------get covergae for intervals which actually have ---
    #--------coverage but no reads covering--------------------
    #only executing when real mate inner dist larger than 0
    interval_l = leftReadsMaxCor[1]
    interval_r = rightReadsMinCor[0]
    overlap = 1 if interval_l > interval_r else 0  
    if extend:
        for pos in xrange(interval_l,interval_r):
            #if pos not in wigDict:
            #    wigDict[pos] = {}
            if xs not in wigDict[pos]:
                wigDict[pos][xs] = [0,-1]
            else:
                wigDict[pos][xs][1] -= 1
            #----------------END adding wigDict-----------------
    #-------END coverage for interval regions---------------
    #if leftReadsMaxCor[1] < rightReadsMinCor[0]:
    #-first get coverage for really mapped regions----
    uniqMappedPosD = {} #urgly unefficient soluable methods
    for reads in pairL:
        for region in reads[2]:
            for pos in xrange(region[0], region[1]):
                if overlap:
                    if pos not in uniqMappedPosD:
                        uniqMappedPosD[pos]=''
                    else:
                        continue #ignore enlarge coverage for position
                                 #sequences twice by one frag
                #if pos not in wigDict:
                #    wigDict[pos] = {}
                if xs not in wigDict[pos]:
                    wigDict[pos][xs] = [1,0]
                else:
                    wigDict[pos][xs][0] += 1
                #------------end adding to wigDict---------
            #----------------end one mapped fragments---
        #----------------end getting each mapped fragments---
    #-----------------end mapped coverage of tow reads--------------------------
    #return wigDict
#----------END computeWigDict---------------

def extendWigDict(wigDict, posKeyL, exonDictChr):
    '''Give coverage to interval regions based on the following two
    conditions.
    1.The interval located at an expressed exons.
    2.The interval has covered by other reads. If this region can not
    be sequenced randomly, it is no need to add coverage for it. 
    Actually no much changes whether add or not.

    Note: This is mostly a rough estimation of coverage especially when
    alternative splicing happened.
    '''
    #xonDictChr = {'+':[[1, 100], [], ...], '-':[]}
    #wigDict = {pos:{+:[+,+_e], '-':[ -,-_e]}}
    positiveExonL = exonDictChr['+']
    positiveExonL.sort()
    if '-' in exonDictChr:
        negativeExonL = exonDictChr['-']
        negativeExonL.sort()
    for pos in posKeyL:
        if wigDict[pos]['+'][1] < 0 and wigDict[pos]['+'][0] > 0:
            for exonR in positiveExonL[:]:
                if exonR[0] <= pos <= exonR[1]:
                    wigDict[pos]['+'][0] += (-1) * wigDict[pos]['+'][1]
                    break
                elif exonR[1] < pos:
                    positiveExonL.remove(exonR) #delete skipped exons
                elif pos < exonR[0]: 
                    break
            #---------test if exon-----------------
        #---------------END positive or unstrand----------
        if '-' in wigDict[pos]:
            if wigDict[pos]['-'][1] < 0 and wigDict[pos]['-'][0] > 0:
                for exonR in negativeExonL[:]:
                    if exonR[0] <= pos <= exonR[1]:
                        wigDict[pos]['-'][0] += (-1) * wigDict[pos]['-'][1]
                        break
                    elif exonR[1] < pos:
                        negativeExonL.remove(exonR) #delete skipped exons
                    elif pos < exonR[0]: #not in Exon 
                        break
                #---------test if exon-----------------
        #---------------END negative or unstrand----------
    #--------------END all position---------------------
#--------NED extendWigDict---------------------

def extendWigDictSE(wigDict):
    global DEBUG
    for pos in wigDict.keys():
        for xs in wigDict[pos].keys():
            real = wigDict[pos][xs][0]
            extend = wigDict[pos][xs][1]
            if DEBUG:
                print >>sys.stderr, "%s\t%s\t%s" % (pos, real, extend)
            if real > 0 and extend < 0:
                wigDict[pos][xs][0] = real - extend

#---------------------------------------------

def outputWigDict(wigDict, chr, posL, lt):
    '''
    Output variableStep wig.
    '''
    #wigDict = wigDict
    if lt == 'fr-unstranded':
        print 'variableStep chrom=%s' % chr
        for i in posL:
            posCov = wigDict[i]['+'][0]
            if posCov:
                print "%d\t%d" % (i,posCov)
    else:
        print "#Columns: Pos, Positive Strand, Negative Strand"
        print 'variableStep chrom=%s' % chr
        for i in posL:
            posCov = wigDict[i]['+'][0]
            negCov = wigDict[i]['-'][0]
            if posCov or negCov:
                print "%d\t%d\t%d" % (i,posCov,negCov)
#--------------output wig by chr----------------------

def outputWigDictTest(wigDict, chr, posL, lt):
    '''
    Output variableStep wig.
    '''
    if lt == 'fr-unstranded':
        print 'variableStep chrom=%s' % chr
        for i in posL:
            posCov = wigDict[i]['+'][0]
            interval = wigDict[i]['+'][1] 
            print "%d\t%d\t%d" % (i,posCov,interval)
    else:
        print "#Columns: Pos, Positive Strand, Negative Strand"
        print 'variableStep chrom=%s' % chr
        for i in posL:
            posCov = wigDict[i]['+'][0]
            negCov = wigDict[i]['-'][0]
            if posCov or negCov:
                print "%d\t%d\t%d" % (i,posCov,negCov)
#--------------output wig by chr----------------------

def outputWigDictSE(wigDict, chr, posL, lt):
    '''
    Output variableStep wig.
    '''
    #wigDict = wigDict
    if lt == 'fr-unstranded':
        print 'variableStep chrom=%s' % chr
        for i in posL:
            posCov = wigDict[i]['+'][0]
            if posCov:
                print "%d\t%d" % (i,posCov)
    else:
        print "#Columns: Pos, Positive Strand, Negative Strand"
        print 'variableStep chrom=%s' % chr
        for i in posL:
            posCov = wigDict[i]['+'][0]
            negCov = wigDict[i]['-'][0]
            if posCov or negCov:
                print "%d\t%d\t%d" % (i,posCov,negCov)
#--------------output wig by chr----------------------




def readExonRegFromGTF(gtf, lt, type='PE'):
    '''
    Get exon regions from GTF for each chromosome. It is such a greedy
    process that if one NT is exon in one transcript, it will be
    considered as exon [This does not always mean it will have
    coverage if interval value less than 0].

    GTF: 1-based numbering, both closed
    '''
    exonDict = {} #{chr:{'+':[[exon_s,exon_e], [], ...], '-':[[],]}} 
                  # 1-based both closed
    for line in open(gtf):
        lineL = line.split('\t', 7)
        if lineL[2] == 'exon':
            xs = lineL[6] if lt != 'fr-unstranded' else '+'
            chr = lineL[0]
            if chr not in exonDict:
                exonDict[chr] = {}
            if xs not in exonDict[chr]:
                exonDict[chr][xs] = []
            exonDict[chr][xs].append([int(lineL[3]), int(lineL[4])])
        #-----------end exon-----------
    #--------End reading file------------
    for chr in exonDict.keys():
        for xs in exonDict[chr].keys():
            exonDict[chr][xs].sort()
    return exonDict
#-----------END read exons from GTF----------

def readExonRegFromGTFSE(gtf, lt):
    '''
    Get exon regions from GTF for each chromosome. It is such a greedy
    process that if one NT is exon in one transcript, it will be
    considered as exon [This does not always mean it will have
    coverage if interval value less than 0].

    GTF: 1-based numbering, both closed
    '''
    exonDict = {} #{chr:{'+':{pos, 1}, '-':[[],]}} 
                  # 1-based both closed
    for line in open(gtf):
        lineL = line.split('\t', 7)
        if lineL[2] == 'exon':
            xs = lineL[6] if lt != 'fr-unstranded' else '+'
            chr = lineL[0]
            if chr not in exonDict:
                exonDict[chr] = {}
            if xs not in exonDict[chr]:
                exonDict[chr][xs] = {}
            for i377 in xrange(int(lineL[3]), int(lineL[4])):
                if i377 not in exonDict[chr][xs]:
                    exonDict[chr][xs][i377]=1
        #-----------end exon-----------
    #--------End reading file------------
    return exonDict
#-----------END read exons from GTF----------



###def computeWigDictSE(wigDict, regionL, start, xs, 
###    extend, flag, exonDictChr=[]):
###    global DEBUG
###    if DEBUG:
###        print >>sys.stderr, "RegionL: ", regionL
###    #------save mapped reads-------------
###    countMappedLen = 0
###    for posL in regionL:
###        for pos in xrange(posL[0], posL[1]):
###            countMappedLen += 1
###            if xs not in wigDict[pos]:
###                wigDict[pos][xs] = [1,0]
###            else:
###                wigDict[pos][xs][0] += 1
###        #----------finish on region-----------
###    #-----------finish all regions-------
###    #if DEBUG:
###    #    outputWigDictTest(wigDict, 'test', posL, lt)
###    if extend and extend > countMappedLen:
###        diff = extend - countMappedLen
###        if DEBUG:
###            print >>sys.stderr, "Diff is %d" % diff
###        #-------------------------------------
###        #exonDict = {} #{chr:{'+':[[exon_s,exon_e], [], ...], '-':[[],]}} 
###                  # 1-based both closed
###        positiveExonL = exonDictChr['+']
###        if '-' in exonDictChr:
###            negativeExonL = exonDictChr['-']
###        #-----------------------------------
###        add = 1
###        findExon = 0
###        if xs == '+':
###            #Origin from positive strand or unstranded reads map to
###            #positive strand.
###            if flag & 16 == 0: 
###                if DEBUG:
###                    print >>sys.stderr, "Flag %d" % flag
###                for exonR in positiveExonL[:]:
###                    if diff == 0:
###                        break
###                    if findExon == 1:
###                        if DEBUG:
###                            print >>sys.stderr, "Current exon:", exonR
###                            print >>sys.stderr, "Next position:", exonR[0]
###                        add = 0
###                        inner_pos = exonR[0]
###                        inner_pos = inner_pos+add
###                        while diff != 0:
###                            if inner_pos > exonR[1]:
###                                break
###                            if xs not in wigDict[inner_pos]:
###                                wigDict[inner_pos][xs] = [0,-1]
###                            else:
###                                wigDict[inner_pos][xs][1] -= 1
###                            inner_pos += 1
###                            diff -= 1
###                            #-------end while-----------------
###                        #-------------------------------------
###                    #------Met the first exon-----------
###                    else:
###                        inner_pos_4 = pos + add
###                        if exonR[0] <= inner_pos_4 <= exonR[1]:
###                            if DEBUG:
###                                print >>sys.stderr, "First Met exon:", exonR
###                                print >>sys.stderr, "Next position:", pos+add
###                            findExon = 1
###                            while inner_pos_4 <= exonR[1]:
###                                if xs not in wigDict[inner_pos_4]:
###                                    wigDict[inner_pos_4][xs] = [0,-1]
###                                else:
###                                    wigDict[inner_pos_4][xs][1] -= 1
###                                #add += 1
###                                inner_pos_4 += 1
###                                diff -= 1
###                                if diff == 0:
###                                    break
###                            #-------end while-----------------
###                        #---do not locate at known exons-----------------
###                        if add == 1 and exonR[0] > pos+add:
###                            if DEBUG:
###                                print >>sys.stderr, "Locate at intron"
###                            break
###                    #---End met the first exon---------------
###                #---------test if exon-----------------
###            #Origin from positive strand or unstranded reads map to
###            #negative strand.
###            else:
###                if DEBUG:
###                    print >>sys.stderr, "Flag %d" % flag
###                for exonR in positiveExonL[::-1]:
###                    if diff == 0:
###                        break
###                    if findExon == 1:
###                        #add = 0
###                        inner_pos = exonR[1]
###                        if DEBUG:
###                            print >>sys.stderr, "Current exon:", exonR
###                            print >>sys.stderr, "Next position:", inner_pos
###                        while diff != 0:
###                            if inner_pos < exonR[0]:
###                                break
###                            if xs not in wigDict[inner_pos]:
###                                wigDict[inner_pos][xs] = [0,-1]
###                            else:
###                                wigDict[inner_pos][xs][1] -= 1
###                            #add -= 1
###                            inner_pos -= 1
###                            diff -= 1
###                            #-------end while-----------------
###                        #-------------------------------------
###                    #---------Try to find the first exon-----
###                    else:
###                        newstart = start - 1
###                        if exonR[0] <= newstart <= exonR[1]:
###                            findExon = 1
###                            if DEBUG:
###                                print >>sys.stderr, "First Met exon:", exonR
###                                print >>sys.stderr, "Next position:", newstart
###                            while newstart >= exonR[0]:
###                                if xs not in wigDict[newstart]:
###                                    wigDict[newstart][xs] = [0,-1]
###                                else:
###                                    wigDict[newstart][xs][1] -= 1
###                                #add += 1
###                                newstart -= 1
###                                diff -= 1
###                                if diff == 0:
###                                    break
###                            #-------end while-----------------
###                        #-------------------------------------
###                        if add == 1 and start-add >exonR[1]:
###                            if DEBUG:
###                                print >>sys.stderr, "Locate at intron"
###                            break
###                    #-------------------------------------
###                #---------test if exon-----------------
###        #---------------END positive or unstrand----------
###        elif xs == '-':
###            if flag & 16 == 0: 
###                for exonR in negativeExonL[:]:
###                    if diff == 0:
###                        break
###                    if findExon == 1:
###                        add = 0
###                        inner_pos = exonR[0] + 1
###                        while diff != 0:
###                            if inner_pos > exonR[1]:
###                                break
###                            if xs not in wigDict[inner_pos]:
###                                wigDict[inner_pos][xs] = [0,-1]
###                            else:
###                                wigDict[inner_pos][xs][1] -= 1
###                            #add += 1
###                            inner_pos += 1
###                            diff -= 1
###                            #-------end while-----------------
###                        #-------------------------------------
###                    else:
###                        inner_pos = pos + 1
###                        if exonR[0] <= inner_pos <= exonR[1]:
###                            findExon = 1
###                            while inner_pos <= exonR[1]:
###                                if xs not in wigDict[inner_pos]:
###                                    wigDict[inner_pos][xs] = [0,-1]
###                                else:
###                                    wigDict[inner_pos][xs][1] -= 1
###                                #add += 1
###                                inner_pos += 1
###                                diff -= 1
###                                if diff == 0:
###                                    break
###                            #-------end while-----------------
###                        #---do not locate at known exons-----------------
###                        if add == 1 and exonR[0] > pos+add:
###                            if DEBUG:
###                                print >>sys.stderr, "Locate at intron"
###                            break
###                    #-------------------------------------
###                #---------test if exon-----------------
###            else:
###                for exonR in negativeExonL[::-1]:
###                    if diff == 0:
###                        break
###                    if findExon == 1:
###                        #add = 0
###                        inner_pos = exonR[1]
###                        while diff != 0:
###                            if inner_pos < exonR[0]:
###                                break
###                            if xs not in wigDict[inner_pos]:
###                                wigDict[inner_pos][xs] = [0,-1]
###                            else:
###                                wigDict[inner_pos][xs][1] -= 1
###                            #add -= 1
###                            inner_pos -= 1
###                            diff -= 1
###                            #-------end while-----------------
###                        #-------------------------------------
###                    else:
###                        inner_pos_5 = start - 1
###                        if exonR[0] <= inner_pos_5 <= exonR[1]:
###                            findExon = 1
###                            while inner_pos_5 >= exonR[0]:
###                                if xs not in wigDict[inner_pos_5]:
###                                    wigDict[inner_pos_5][xs] = [0,-1]
###                                else:
###                                    wigDict[inner_pos_5][xs][1] -= 1
###                                #add += 1
###                                inner_pos_5 -= 1
###                                diff -= 1
###                                if diff == 0:
###                                    break
###                            #-------end while-----------------
###                        #-----locate at intron-------------
###                        if add == 1 and start-add >exonR[1]:
###                            if DEBUG:
###                                print >>sys.stderr, "Locate at intron"
###                            break
###                    #---------find the first exon-----------
###                #---------test if exon-----------------
###            #---------------------------
###        #---------------END positive or unstrand----------
###    #--------------END all position---------------------
####---------END SE-----------------------------------


def computeWigDictSE(wigDict, regionL, start, xs, 
    extend, flag, exonDictChr={}):
    global DEBUG
    if DEBUG:
        print >>sys.stderr, "RegionL: ", regionL
    #------save mapped reads-------------
    countMappedLen = 0
    for posL in regionL:
        for pos in xrange(posL[0], posL[1]):
            countMappedLen += 1
            if xs not in wigDict[pos]:
                wigDict[pos][xs] = [1,0]
            else:
                wigDict[pos][xs][0] += 1
        #----------finish on region-----------
    #-----------finish all regions-------
    #if DEBUG:
    #    outputWigDictTest(wigDict, 'test', posL, lt)
    if extend and extend > countMappedLen:
        diff = extend - countMappedLen
        if DEBUG:
            print >>sys.stderr, "Diff is %d" % diff
        #-------------------------------------
        #exonDict = {} #{chr:{'+':[[exon_s,exon_e], [], ...], '-':[[],]}} 
                  # 1-based both closed
        positiveExonD = exonDictChr['+']
        if '-' in exonDictChr:
            negativeExonD = exonDictChr['-']
        #-----------------------------------
        add = 1
        findExon = 0
        if xs == '+':
            #Origin from positive strand or unstranded reads map to
            #positive strand.
            if flag & 16 == 0: 
                if DEBUG:
                    print >>sys.stderr, "Flag %d" % flag
                inner_pos = pos
                if inner_pos in positiveExonD:
                    while diff != 0:
                        inner_pos += 1
                        if inner_pos in positiveExonD:
                            if xs not in wigDict[inner_pos]:
                                wigDict[inner_pos][xs] = [0,-1]
                            else:
                                wigDict[inner_pos][xs][1] -= 1
                            diff -= 1
                        #---------test if exon-----------------
                        if inner_pos == 1000000000:
                            break
                    #---------------------------------------------
                #-------Extend only when reads located at exons----
            #Origin from positive strand or unstranded reads map to
            #negative strand.
            else:
                if DEBUG:
                    print >>sys.stderr, "Flag %d" % flag
                inner_pos = start
                if inner_pos in positiveExonD:
                    while diff != 0:
                        inner_pos -= 1
                        if inner_pos in positiveExonD:
                            if xs not in wigDict[inner_pos]:
                                wigDict[inner_pos][xs] = [0,-1]
                            else:
                                wigDict[inner_pos][xs][1] -= 1
                            diff -= 1
                        #---------test if exon-----------------
                        if inner_pos == -1:
                            break
                    #------------------------------------
                #-------Extend only when reads located at exons----
        #---------------END positive or unstrand----------
        elif xs == '-':
            if flag & 16 == 0: 
                inner_pos = pos
                if inner_pos in negativeExonD:
                    while diff != 0:
                        inner_pos += 1
                        if inner_pos in nagativeExonD:
                            if xs not in wigDict[inner_pos]:
                                wigDict[inner_pos][xs] = [0,-1]
                            else:
                                wigDict[inner_pos][xs][1] -= 1
                            diff -= 1
                        #---------test if exon-----------------
                        if inner_pos == 1000000000:
                            break
                    #---------test if exon-----------------
                #-------Extend only when reads located at exons----
            else:
                inner_pos = start
                if inner_pos in negativeExonD:
                    while diff != 0:
                        inner_pos -= 1
                        if inner_pos in negativeExonD:
                            if xs not in wigDict[inner_pos]:
                                wigDict[inner_pos][xs] = [0,-1]
                            else:
                                wigDict[inner_pos][xs][1] -= 1
                            diff -= 1
                        #---------test if exon-----------------
                        if inner_pos == -1:
                            break
                    #------------------------------------
                #-------Extend only when reads located at exons----
            #---------------------------
        #---------------END nagative----------
    #--------------END all position---------------------
#---------END SE-----------------------------------




def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    cigarP = re.compile('([0-9]+)([A-Z])')
    file = options.filein
    readsType = options.seq_Type
    lt= options.lt
    extend = options.extend
    gtf = options.gtf
    nt = options.nt
    cs = options.chromSize
    verbose = options.verbose
    global DEBUG
    DEBUG = options.debug
    #+_p means number of directly mapped reads in positive strand of this postion
    #+_e means number of potential interval reads mapped to positive
    #strand of thi spostion
    #- means negative strand.
    #wigDict = {} #dict = {pos:{+_p:[+,+_e], '-':[-_p,-_e]}}
    wigDict = collections.defaultdict(dict)
    wigDictSE = {}
    pairDict = {}
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    #-------------open GTF-----------
    if nt == 'RNA' and extend:
        if verbose:
            print >>sys.stderr, "--Begin reading GTF---%s" \
                    % strftime(timeformat, localtime())
        if gtf:
            if readsType == 'PE':
                exonDict = readExonRegFromGTF(gtf,lt)   
            elif readsType == 'SE':
                exonDict = readExonRegFromGTFSE(gtf, lt)
        else:
            print >>sys.stderr, "GTF file is needed."
            sys.exit(1)
        if verbose:
            print >>sys.stderr, "--Finish reading GTF---%s" \
                    % strftime(timeformat, localtime())
    #---------------------------------
    chr = ''
    count = 0
    #read in sam file
    for line in fh:
        if verbose:
            count += 1
            if count % 100000 == 0:
                print >>sys.stderr, "------Finish %d reads----%s" \
                    % (count, strftime(timeformat, localtime()))
        #---------------------------------------------------
        lineL = line.strip().split("\t")
        name = lineL[0]
        if DEBUG:
            print >>sys.stderr, "---Processing %s" % name
        flag = int(lineL[1])
        #----------------------START output one chr----------------
        if chr and chr != lineL[2]:
            posL = wigDict.keys()
            posL.sort()
            if readsType == 'PE' and extend and chr in exonDict:
                if verbose:
                    print >>sys.stderr,\
                        "--Begin extend PE reads for %s--%s" \
                        % (chr, strftime(timeformat, localtime()))
                extendWigDict(wigDict, posL, exonDict[chr])
                exonDict.pop(chr)
            #--Extend single end sequencing of RNA-Seq 
            elif readsType == 'SE' and extend:
                if verbose:
                    print >>sys.stderr,\
                        "--Begin extend SE reads for %s--%s" \
                        % (chr, strftime(timeformat, localtime()))
                extendWigDictSE(wigDict)
            if verbose:
                print >>sys.stderr,\
                    "--Begin output wig for %s--%s" \
                        % (chr, strftime(timeformat, localtime()))
            if DEBUG:
                outputWigDictTest(wigDict, chr, posL, lt)
            else:
                outputWigDict(wigDict, chr, posL, lt)
            wigDict = {}
            wigDict = collections.defaultdict(dict)
        #----------------------END output one chr----------------
        chr = lineL[2]
        start = int(lineL[3]) ##sam and wig are 1-based
        cigar = lineL[5] 
        regionL = computeRegion(start,cigarP.findall(cigar),name) 
        if lt != 'fr-unstranded':
            xs = [i[-1] for i in lineL[11:] if i.startswith('XS:A:')]
        else:
            xs = '+'
        if readsType == 'SE':
            if extend:
                computeWigDictSE(wigDict, regionL, start, xs, \
                    extend, flag, exonDict[chr])
            else:
                for posL in regionL:
                    for pos in xrange(posL[0], posL[1]):
                        if xs not in wigDict[pos]:
                            wigDict[pos][xs] = [1,0]
                        else:
                            wigDict[pos][xs][0] += 1
                    #----------finish on region-----------
                #-----------finish all regions-------
            #---------Extend or not----------------------
        #---------------END SE----------------------
        elif readsType == 'PE':
            if flag & 0x2 == 2 and lineL[6] == '=': #properly paired
                if name not in pairDict:
                    pairDict[name] = [[chr,flag,regionL,xs]]
                else:
                    pairDict[name].append([chr,flag,regionL,xs])
                    computeWigDict(wigDict, pairDict[name], extend)
                    pairDict.pop(name) #here only pop is right. clear
                    #out whole dict will give wrong results when two
                    #paired reads are nor neighbors.
                #------------------------------
            elif flag & 0x2 == 0 or lineL[6] != '=': #unproperly paired
                for posL in regionL:
                    for pos in xrange(posL[0], posL[1]):
                        #if pos not in wigDict:
                        #    wigDict[pos] = {}
                        if xs not in wigDict[pos]:
                            wigDict[pos][xs] = [1,0]
                        else:
                            wigDict[pos][xs][0] += 1
                    #----------finish on region-----------
                #-----------finish all regions-------
            #--------END unproperly paired---------
            else:
                print >>sys.stderr, "Unknown %s" % line
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    #----last chromosome-------------------------------
    if wigDict:
        posL = wigDict.keys()
        posL.sort()
        if readsType == 'PE' and extend and chr in exonDict:
            if verbose:
                print >>sys.stderr,\
                    "--Begin extend PE reads for %s--%s" \
                    % (chr, strftime(timeformat, localtime()))
            extendWigDict(wigDict, posL, exonDict[chr])
        elif readsType == 'SE' and extend:
            if verbose:
                print >>sys.stderr,\
                    "--Begin extend SE reads for %s--%s" \
                    % (chr, strftime(timeformat, localtime()))
            extendWigDictSE(wigDict)
        if verbose:
            print >>sys.stderr,\
                "--Begin output wig for %s--%s" \
                    % (chr, strftime(timeformat, localtime()))
            #exonDict.pop(chr)
        if DEBUG:
            outputWigDictTest(wigDict, chr, posL, lt)
        else:
            outputWigDict(wigDict, chr, posL, lt)

    #--------------------------------------------------
    if verbose:
        print >>sys.stderr,\
            "--Successful %s--%s"
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


