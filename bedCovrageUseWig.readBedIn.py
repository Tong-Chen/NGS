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
Compute the average read coverage of given bed file using wig. 

Bed file (other lines are not allowed):
chr7    52823164        52823749        0610005C13Rik   -
chr7    52826355        52826562        0610005C13Rik   -
chr7    52829782        52829892        0610005C13Rik   -
chr7    52829977        52830147        0610005C13Rik   -
chr7    52830496        52830546        0610005C13Rik   -
chr5    31351012        31351129        0610007C21Rik   +
chr5    31351834        31351953        0610007C21Rik   +
chr5    31354569        31354641        0610007C21Rik   +
chr5    31354834        31354906        0610007C21Rik   +
chr5    31355135        31355257        0610007C21Rik   +
chr5    31356333        31356431        0610007C21Rik   +

The forth column may have same name or different names. Lines
with same names will be merged together to compute sum,mean,median,
max,min. Howeveer, Only the last line will be outputed.

Wig file (all regions belong to one chromosome must be under one
'variableStep' type line):
track type=wig nanme=""
variableStep chrom=chr1 span=10
3001321 0.023209
3001331 0.023209
3001341 0.023209
3001351 0.023209
3001361 0.023209
3001371 0.023209
3001381 0.023209
""""
'''

import collections
from numpy import mean,median,max,min,sum,argsort
from numpy import array as np_array
from numpy import append as np_append
from array import array
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
    usages = "%prog -i bed -w wig -o operator -s True"
    parser = OP(usage=usages)
    parser.add_option("-i", "--bed", dest="bed",
        metavar="BED_REGION", help="Regions in bed file format, \
Bed will be read in memory wholely. - can be used as STDIN.\
When -m is TRUE, no duplicate names allowed in bed. Otherwise \
neighbor lines \
with same name in the forth column will be taken as one regions.***")
    parser.add_option("-w", "--wig", dest="wig",
        metavar="WIG", help="Regions in wig file format. Each position \
of same chromsome in wig must be sorted numerically. All legal wig \
obeys this rule. Tested for variable step and fixed step wig. \
All positions of one chromosome must under one header line when using \
variable step wig. For fixed step wig, all separated regions of one \
chromosome must near each other and should not be interrupted by \
regions in other chromosomes.***")
    parser.add_option("-o", "--op", dest="op",
        metavar="OPERATOR", help="Several choice, sum,mean,median,max,min.\
Multiple ones can be given in ',' connected in formats like \
<mean, max>. This parameter is exclusive with -m.")
    parser.add_option("-m", "--maximum-pos", dest="mp",
        metavar="1/0", default=0, help="Get the position with the largest \
value. If multiple maximum positions are found, the one nearest the \
middle will be chosed. This parameter is exclusive with -o. Default 0 \
means do not get max value.")
    parser.add_option("-s", "--strand", dest="strand",
        metavar="1/0", default=0, help="When 1 is given, compute \
strand specific coverage for bed regions. This assumes, the seond \
column in wig is positive strand while the third column in wig is \
nagative strand. Default 0, means unstranded." )
    parser.add_option("-n", "--name", dest="name", default=1,
        metavar="All column/Name column", help="Default the \
column output before coverage data is the forth column of bed file. If\
position information or all columns are also required,  please give <0>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.bed != None, "Region file in bed format is needed -i"
    assert options.wig != None, "Coverage file in wig format is needed -i"
    return (options, args)
#--------------------------------------------------------------------


def readBed(bed_fh):
    '''
    Read bed into a dicat indexed by chromosome name.

    additional = [0, 0, 0] length of additional is the same as length for
    opL.

    bedDict = {chr:[[chr, int_start, int_end, name, other, additional], \
        [], []]}
    '''
    bedDict = {}
    for line in bed_fh:
        lineL = line.strip().split('\t')
        #lineL[1] = int(lineL[1])
        #lineL[2] = int(lineL[2])
        chr = lineL[0]
        #lineL.extend(additional)
        if chr not in bedDict:
            bedDict[chr] = []
        bedDict[chr].append(lineL)
    #-------------------------------
    chrRegion= {} # chrRegion = {chr:(start, end)}
    for chr, valueL in bedDict.items():
        start = [int(i[1]) for i in valueL]
        end   = [int(i[2]) for i in valueL]
        chrRegion[chr] = (min(start), max(end))
    return bedDict, chrRegion
#-------------------------------------------------

#--------------------------------------------
def computeCoverage(bedDict, wig, opL, name_mode, strand, chrRegion, mp):
    #array_i = array
    chr = ''
    span = 0
    step = 0
    #wigDict = collections.defaultdict(dict)
    #wigDict = {} #unstrand wigDict = {pos:value}
                #strand wigDicr = {pos:[pos value, neg value]}
    bedDictChr = []
    valDict = {}
    pos_fixed = 0
    overlappedChr = []
    for i in open(wig):
        if i.startswith('track'):
            continue
        elif i.startswith('#'):
            continue
        elif i.startswith('browse'):
            continue
        elif i.startswith('variableStep'):
            saveWig = 1
            if chr and chr in bedDict:
                outputCoverage(wigChrDict,bedDict[chr],opL,name_mode,\
                    strand, mp)
                bedDict.pop(chr)
            wigChrDict = {}
            pos_fixed = 0
            chromi = i.rfind('chrom=')
            assert chromi != -1, "Wrong format no chr %s" % i
            chr = i[chromi+6:].strip().split()[0]
            overlappedChr.append(chr)
            if chr not in bedDict:
                saveWig = 0
            else:
                spani = i.rfind("span=")
                span = 1 if spani == -1 else \
                    int(i[spani+5:].strip().split()[0])
                pos = 1
                neg = 2
                min_start = chrRegion[chr][0]
                max_end   = chrRegion[chr][1]
        #-unckecked for less of data-------------
        elif i.startswith("fixedStep"):
            saveWig = 1
            #if chr and chr in bedDict and wigChrDict:
            chromi = i.rfind('chrom=')
            assert chromi != -1, "Wrong format no chr %s" % i
            newchr = i[chromi+6:].strip().split()[0]
            if chr and chr != newchr and chr in bedDict:
                outputCoverage(wigChrDict,bedDict[chr],opL,\
                    name_mode,strand,mp)
                bedDict.pop(chr)
                wigChrDict = {}
                overlappedChr.append(chr)
            if not chr:
                wigChrDict = {}
            chr = newchr
            if chr not in bedDict:
                saveWig = 0
            else:
                starti = i.rfind('start=')
                assert starti != -1, "Must have start %s" % i
                pos_fixed = int(i[starti+6:].strip().split()[0])
                assert pos_fixed > 0, \
                    "fixedStep must have start bigger than 0 %s" % i
                stepi = i.rfind('step=')
                assert stepi != -1, "Must have step %s" % i
                step = int(i[stepi+5:].strip().split()[0])
                start = pos_fixed - step - 1 #feasible add <step> each
                spani = i.rfind('span=')
                span = 1 if spani == -1 else \
                    int(i[spani+5:].strip().split()[0])
                pos = 0
                neg = 1
                min_start = chrRegion[chr][0]
                max_end   = chrRegion[chr][1]
            #--------------------------------------------------
        else:
            if not saveWig:
                continue
            lineL = i.strip().split()
            if pos_fixed: #fixedStep each position is the last
                        #position plus step. 
                start += step
                end = start + span   
            else:
                start = int(lineL[0])-1
                end = start + span
            #-------------------------------------------
            #skip coverage for regions not located in bed
            if end <= min_start:
                continue
            if start > max_end:
                continue
            #---------------------------------------------
            if strand:
                for position in xrange(start,end):
                    wigChrDict[position] = \
                        [float(lineL[pos]),float(lineL[neg])] 
            else:
                for position in xrange(start,end):
                    wigChrDict[position] = float(lineL[pos]) 
            #-------------------------------------------------
        #--------end processing one line ----------------------
    #----------------END reading whole file-------------------
    if chr and chr in bedDict:
        outputCoverage(wigChrDict,bedDict[chr],opL,\
            name_mode,strand,mp)
        bedDict.pop(chr)
    return overlappedChr
#---------------------------------------------------

def outputCoverage(wigChrDict,bedLineL,opL,name_mode,strand,mp):
    label = ''
    valueL = np_array([])
    other = ''
    for lineL in bedLineL:
        if label and label != lineL[3]:
            #---------------get the position with maximum value-------
            if mp:
                sort_indexL = argsort(valueL)
                #print valueL
                #print sort_indexL
                maxP = sort_indexL[-1]
                #print 'ORiginal maxP:', maxP
                diff = abs(sort_indexL[-1]-mid)
                maxV = valueL[maxP]
                #print 'Total length:', length
                for i in xrange(-2,length,-1):
                    if valueL[sort_indexL[i]] == maxV:
                        if diff > abs(sort_indexL[i]-mid):
                            diff = abs(sort_indexL[i]-mid)
                            maxP = sort_indexL[i]
                            maxV = valueL[maxP]
                            #print 'Iterated maxP:', maxP
                    else:
                        break
                #----------------------------
                if other:
                    print '%s\t%s\t%s\t%s\t%s\t%s' % (chr, str(maxP+start),
                            str(maxP+start+1), name, str(maxV),
                            '\t'.join(other))
                    other = ''
                else:
                    print '%s\t%s\t%s\t%s\t%s' % (chr, str(maxP+start),
                            str(maxP+start+1), name, str(maxV))
            #---------------get the position with maximum value-------
            else:
                print "%s\t%s" % (name, \
                    '\t'.join([str(op(valueL)) for op in opL]))
            valueL = np_array([])
        #----------output ----------------------------------
        #---------------begin process--------------------
        chr = lineL[0]
        start = int(lineL[1])
        end   = int(lineL[2])
        length = start - end - 1 #used for xrange, -1 means include
                                #start-end
        mid = (end-start)/2.0
        label = lineL[3]
        if len(lineL) >= 6:
            other = '\t'.join(lineL[5:])
        if strand:
            strand_in = lineL[5]
            strand_num = 0 if strand_in == '+' else 1
            if name_mode:
                name = ''.join(label, '@', strand_in)
            else:
                name = '\t'.join(lineL)
            valueL = np_append(valueL, \
                [wigChrDict.get(i,[0,0])[strand_num] for i in xrange(start, end)])
        else:
            if name_mode:
                name = label
            else:
                name = '\t'.join(lineL)
            valueL = np_append(valueL, \
                [wigChrDict.get(i,0) for i in xrange(start, end)])
    #---------------------------------------------------------------------
    if label:
        #---------------get the position with maximum value-------
        if mp:
            sort_indexL = argsort(valueL)
            #print valueL
            #print sort_indexL
            maxP = sort_indexL[-1]
            #print 'ORiginal maxP:', maxP
            diff = abs(sort_indexL[-1]-mid)
            maxV = valueL[maxP]
            #print 'Total length:', length
            for i in xrange(-2,length,-1):
                if valueL[sort_indexL[i]] == maxV:
                    if diff > abs(sort_indexL[i]-mid):
                        diff = abs(sort_indexL[i]-mid)
                        maxP = sort_indexL[i]
                        #print maxP
                else:
                    break
            #----------------------------
            if other:
                print '%s\t%s\t%s\t%s\t%s\t%s' % (chr, str(maxP+start),
                        str(maxP+start+1), name, str(maxV),
                        '\t'.join(other))
            else:
                print '%s\t%s\t%s\t%s\t%s' % (chr, str(maxP+start),
                        str(maxP+start+1), name, str(maxV))
            #print "%s\t%s" % (name, str(maxP+start))
        #---------------get the position with maximum value-------
        else:
            print "%s\t%s" % (name, \
                '\t'.join([str(op(valueL)) for op in opL]))
        valueL = np_array([])

#---------------------------------------------------
#
#def output_new(bedLineL, valList, opL, opDict, name_mode, strand):
#    regionL = bedLineL[2]-bedLineL[1]
#    bedLineL[1] = str(bedLineL[1])
#    bedLineL[2] = str(bedLineL[2])
#    #print bedLineL
#    name = bedLineL[3]
#    LenvalList = len(valList)
#    diff = regionL - LenvalList
#    if strand:
#        valList.extend([(0,0) for j in range(diff)])
#        posValL = [vp[0] for vp in valList]
#        negValL = [vp[1] for vp in valList]
#        output = '\t'.join([\
#            '\t'.join([str(opDict[op](posValL)) for op in opL]), \
#            '\t'.join([str(opDict[op](negValL)) for op in opL])])
#    else:
#        valList.extend([0 for j in range(diff)])
#        output = '\t'.join([str(opDict[op](valList)) for op in opL])
#    #----------------------------------------------
#    if name_mode:
#        print "%s\t%s" % (name, output)
#    else:
#        print "%s\t%s" % ('\t'.join(bedLineL), output)
#
##--------------------------------------------------------------------
#def output(bedDictChr, valDict, opL, opDict, name_mode, strand):
#    for bedLineL in bedDictChr:
#        regionL = bedLineL[2]-bedLineL[1]
#        bedLineL[1] = str(bedLineL[1])
#        bedLineL[2] = str(bedLineL[2])
#        #print bedLineL
#        name = bedLineL[3]
#        valList = valDict[name]
#        LenvalList = len(valList)
#        diff = regionL - LenvalList
#        if strand:
#            valList.extend([(0,0) for j in range(diff)])
#            posValL = [vp[0] for vp in valList]
#            negValL = [vp[1] for vp in valList]
#            output = '\t'.join([\
#                '\t'.join([str(opDict[op](posValL)) for op in opL]), \
#                '\t'.join([str(opDict[op](negValL)) for op in opL])])
#        else:
#            valList.extend([0 for j in range(diff)])
#            output = '\t'.join([str(opDict[op](valList)) for op in opL])
#        #----------------------------------------------
#        if name_mode:
#            print "%s\t%s" % (name, output)
#        else:
#            print "%s\t%s" % ('\t'.join(bedLineL), output)
#    #--------------------------------------------------------------
#
#def outputWig(wigDict, strand):
#    keyL = wigDict.keys()
#    keyL.sort()
#    for key in keyL:
#        chrD = wigDict[key]
#        chrD_keys = chrD.keys()
#        chrD_keys.sort()
#        print key
#        for i in chrD_keys:
#            print "%d\t%f" % (i,chrD[i])
##----------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    bed = options.bed
    wig = options.wig
    strand = options.strand
    verbose = options.verbose
    debug = options.debug
    mp = options.mp
    #wigDict = readWig(wig, strand)
    #addL = [0 for i in opL]
    name_mode = int(options.name)
    #-----------------------------------
    opDict = {'mean':mean, 'median':median, \
            'max':max, 'min':min, 'sum':sum}
    if options.op:
        opL = [opDict[i] for i in options.op.split(',')]
    else:
        opL = []
        assert options.mp
    if bed == '-':
        fh = sys.stdin
    else:
        fh = open(bed)
    #--------------------------------
    if options.op:
        if name_mode:
            print "#name\t%s" % '\t'.join([i for i in options.op.split(',')])
        else:
            print "#%s" % '\t'.join([i for i in options.op.split(',')])
    bedDict, chrRegion = readBed(fh)
    overlappedChr = computeCoverage(bedDict, wig, opL, name_mode,
            strand, chrRegion, mp)
    #-------output regions which do not exist in wig---
    if options.op:
        valuL = '\t'.join(['0' for i in opL])
        for key in bedDict.keys():
            if key not in overlappedChr:
                for lineL in bedDict[key]:
                    label = lineL[3]
                    if name_mode:
                        if strand:
                            strand_in = lineL[5]
                            name = ''.join(label, '@', strand_in)
                        else:
                            name = label
                    else:
                        name = '\t'.join(lineL)
                    print '%s\t%s' % (name, valuL)
                #--------------------------------------
            #----------------------------------------
        #-----------------------------------------------
    #----close file handle for files-----
    if bed != '-':
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



