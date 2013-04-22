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
from numpy import mean,median,max,min,sum
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
Bed will be read in memory wholely. - can be used as STDIN.***")
    parser.add_option("-w", "--wig", dest="wig",
        metavar="WIG", help="Regions in wig file format. Each position \
of same chromsome in wig must be sorted numerically. All legal wig \
obeys this rule.***")
    parser.add_option("-o", "--op", dest="op",
        metavar="OPERATOR", help="Several choice, sum,mean,median,max,min.\
        Multiple ones can be given in ',' connected formatsi, lke \
        <mean, max>.")
    parser.add_option("-s", "--strand", dest="strand",
        metavar="1/0", default=0, help="When 1 is given, compute \
strand specific coverage for bed regions. This assumes, the seond \
column in wig is positive strand while the third column in wig is \
nagative strand." )
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
    return bedDict
#-------------------------------------------------

#--------------------------------------------
def computeCoverage(bedDict, wig, opL, name_mode, strand):
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
            if chr and chr in bedDict and wigChrDict:
                outputCoverage(wigChrDict,bedDict[chr],opL,name_mode,strand)
                bedDict.pop(chr)
            wigChrDict = {}
            pos_fixed = 0
            chromi = i.rfind('chrom=')
            assert chromi != -1, "Wrong format no chr %s" % i
            chr = i[chromi+6:].strip().split()[0]
            overlappedChr.append(chr)
            if chr not in bedDict:
                saveWig = 0
            spani = i.rfind("span=")
            span = 1 if spani == -1 else \
                int(i[spani+5:].strip().split()[0])
            pos = 1
            neg = 2
        #-unckecked for less of data-------------
        elif i.startswith("fixedStep"):
            saveWig = 1
            if chr and chr in bedDict and wigChrDict:
                outputCoverage(wigChrDict,bedDict[chr],opL,name_mode,strand)
                bedDict.pop(chr)
            wigChrDict = {}
            chromi = i.rfind('chrom=')
            assert chromi != -1, "Wrong format no chr %s" % i
            chr = i[chromi+6:].strip().split()[0]
            overlappedChr.append(chr)
            if chr not in bedDict:
                saveWig = 0
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
    if chr and chr in bedDict and wigChrDict:
        outputCoverage(wigChrDict,bedDict[chr],opL,name_mode,strand)
        bedDict.pop(chr)
    return overlappedChr
#---------------------------------------------------

def outputCoverage(wigChrDict,bedLineL,opL,name_mode,strand):
    label = ''
    valueL = np_array([])
    for lineL in bedLineL:
        start = int(lineL[1])
        end   = int(lineL[2])
        if label and label != lineL[3]:
            print "%s\t%s" % (name, \
                '\t'.join([str(op(valueL)) for op in opL]))
            valueL = np_array([])
        label = lineL[3]
        if strand:
            strand_in = lineL[5]
            strand_num = 0 if strand_in == '+' else 1
            if name_mode:
                name = ''.join(label, '@', strand_in)
            else:
                name = '\t'.join(lineL)
            valueL = np_appedn(valueL, \
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
    #wigDict = readWig(wig, strand)
    #addL = [0 for i in opL]
    name_mode = int(options.name)
    #-----------------------------------
    opDict = {'mean':mean, 'median':median, \
            'max':max, 'min':min, 'sum':sum}
    opL = [opDict[i] for i in options.op.split(',')]
    if bed == '-':
        fh = sys.stdin
    else:
        fh = open(bed)
    #--------------------------------
    if name_mode:
        print "#name\t%s" % '\t'.join([i for i in options.op.split(',')])
    else:
        print "#%s" % '\t'.join([i for i in options.op.split(',')])
    bedDict = readBed(fh)
    overlappedChr = computeCoverage(bedDict, wig, opL, name_mode, strand)
    #-------output regions which do not exist in wig---
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



