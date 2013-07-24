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
Merge two wig files together by supplying parameters like substract,
divide, average or sum.

Bed file:
track line [optional, will be skipped]
#annotation line [optional, will be skipped]
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

The forth column may have same name or different names. Neighbor lines
with same names will be merged together to compute sum,mean,median,
max,min
'''

import collections
import cPickle
from numpy import mean,median,max,min,sum
from numpy import array as np_array
from numpy import append as np_append
from array import array
import sys
import os
from time import localtime, strftime, time
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
verbose = 0

prefix=str(time())

def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    desc = "Output to file"
    usages = "%prog -i bed -w wig -o operator -s True"
    parser = OP(usage=usages)
    parser.add_option("-w", "--wig1", dest="wig1",
        metavar="WIG1", help="A normalized wig file. Each value \
represents the RPM for that point. \
For variable step, one chromosome only allowed to be appear once. \
For fixed step, if one chromosome appears several times in a wig file, they \
must appear continuously.")
    parser.add_option("-W", "--wig2", dest="wig2",
        metavar="WIG2", help="A normalized wig file. Each value \
represents the RPM for that point. \
For variable step, one chromosome only allowed to be appear once. \
For fixed step, if one chromosome appears several times in a wig file, they \
must appear continuously.")
    parser.add_option("-o", "--op", dest="op",
        metavar="OPERATOR", default="substract", help="Several choice, \
divide, substract, average, sum. Multiple ones can be given \
in ',' connected formatsi, like  <divide,substract> \
(without angle brackets). When <divide> and <substract> is given, \
it usually means the first wig divides or substracts the second wig. \
Default 'substract'. ")
    parser.add_option("-s", "--strand", dest="strand",
        metavar="1/0", default=0, help="When 1 is given, compute \
strand specific coverage for bed regions. This assumes, the seond \
column in wig is positive strand while the third column in wig is \
nagative strand." )
    parser.add_option("-p", "--output-prefix", dest="outP",
        help="The prefix for output files")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.wig1 != None, "Wig1 is needed"
    assert options.wig2 != None, "Wig2 is needed"
    return (options, args)
#--------------------------------------------------------------------

def readWig(wig, strand, outputSingle, outputChrList, processL, fileD):
    array_i = array
    chr = ''
    span = 0
    step = 0
    wigDict = collections.defaultdict(dict)
    #wigDict = {} #unstrand wigDict = {pos:value}
                #strand wigDicr = {pos:[pos value, neg value]}
    pos_fixed = 0
    for i in open(wig):
        if i.startswith('track'):
            continue
        elif i.startswith('#'):
            continue
        elif i.startswith('browse'):
            continue
        elif i.startswith('variableStep'):
            if chr and outputSingle:
                saveWigByChr(wigDict, chr, outputChrList)
                wigDict = collections.defaultdict(dict)
            #-------------------------------------------------
            if chr and processL:
                processWig(wigDict, chr, outputChrList, processL,
                    fileD)
                wigDict = collections.defaultdict(dict)
            #-------------------------------------------------------
            pos_fixed = 0
            chromi = i.rfind('chrom=')
            assert chromi != -1, "Wrong format no chr %s" % i
            chr = i[chromi+6:].strip().split()[0]
            spani = i.rfind("span=")
            span = 1 if spani == -1 else \
                int(i[spani+5:].strip().split()[0])
            pos = 1
            neg = 2
        #-unckecked for less of data-------------
        elif i.startswith("fixedStep"):
            chromi = i.rfind('chrom=')
            assert chromi != -1, "Wrong format no chr %s" % i
            newchr = i[chromi+6:].strip().split()[0]
            if chr and chr != newchr and outputSingle:
                saveWigByChr(wigDict, chr, outputChrList)
                wigDict = collections.defaultdict(dict)
            #-------------------------------------------------
            if chr and chr != newchr and processL:
                processWig(wigDict, chr, outputChrList, processL,
                    fileD)
                wigDict = collections.defaultdict(dict)
            #-------------------------------------------------------
            chr = newchr
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
            lineL = i.strip().split()
            if pos_fixed: #fixedStep each position is the last
                        #position plus step. 
                start += step
                end = start + span   
            else:
                start = int(lineL[0])-1
                end = start + span
            #-------------------------------------------
            #if span < 2 and not strand:
            #    wigDict[chr][(start,end)] = float(lineL[pos])
            #elif span < 2 and strand:
            #    wigDict[chr][(start,end)] = array('f', \
            #        [float(lineL[pos]), float(lineL[neg])])
            #else:
            for position in xrange(start, end):
                if not strand:
                    wigDict[chr][position] = float(lineL[pos])
                else:
                    wigDict[chr][position] = array('f', \
                        [float(lineL[pos]), float(lineL[neg])])
                    #--------------------------------------
                #-----------------------------------------
            #-------------------------------------------------
        #--------end processing one line ----------------------
    #----------------END reading whole file-------------------
    #return wigDict
    if chr and outputSingle:
        saveWigByChr(wigDict, chr, outputChrList)
        wigDict = collections.defaultdict(dict)
    #-------------------------------------------------
    if chr and processL:
        processWig(wigDict, chr, outputChrList, processL, fileD)
        wigDict = collections.defaultdict(dict)
    #-------------------------------------------------------
#---------------------------------------------------

def outputWig(wigDict, strand):
    keyL = wigDict.keys()
    keyL.sort()
    for key in keyL:
        chrD = wigDict[key]
        chrD_keys = chrD.keys()
        chrD_keys.sort()
        print key
        for i in chrD_keys:
            print "%d\t%f" % (i,chrD[i])
#----------------------------------------

def saveWigByChr(wigDict, chr, outputChrList):
    file= chr + '.pkl___ct' + prefix
    if file in outputChrList:
        print >>sys.stderr, "Duplicate %s for ourput, \
are you sure all the parts of one chr are together" % file
        sys.exit(1)
    outputChrList.add(chr)
    output = open(file, 'wb')
    cPickle.dump(wigDict[chr], output, 
        protocol=cPickle.HIGHEST_PROTOCOL)
    output.close()


#---------------------------------------------

def processWig(wigDict, chr, outputChrList, processL, fileD):
    '''
    wigDict = {chr:{pos:value}} #here pos is bed format, start with 0
    w1D = {pos:value}
    w2D = {pos:value}
    '''
    #assert chr in outputChrList, \
    #    "Unmatched %s. This may caused by two reasons, first \
    #    you do not have same chromosome list in your \
    #    two wig files. This happens rarely. Second,  \
    #    different parts of your chromosome did not appear \
    #    together." % file
    if chr in outputChrList:
        file = chr + '.pkl___ct' + prefix
        input = open(file, 'rb')
        w1D = cPickle.load(input)
        input.close()
    else:
        print >>sys.stderr, "Unmatched chr %s " % chr
        w1D = {}
        w1D[0] = 0
    #---------------------------
    w2D = wigDict[chr]
    if verbose:
        print >>sys.stderr, "-w wig1"
        print >>sys.stderr, w1D
        print >>sys.stderr, "-W wig2"
        print >>sys.stderr, w2D
    outputChrList.discard(chr)
    os.remove(file)
    w1DKeys = set(w1D.keys())
    w2DKeys = set(w2D.keys())
    posL = list(w1DKeys.union(w2DKeys))
    posL.sort()
    for fh in fileD.values():
        print >>fh, "variableStep chrom=%s" % chr
    for i in posL:
        w1_value = w1D.get(i, 0)
        w2_value = w2D.get(i, 0)
        for op in processL:
            if op == 'substract':
                i_value = w1_value - w2_value
                print >>fileD[op], "%d\t%f" % (i+1, i_value)
            elif op == 'divide':
                if w2_value == 0:
                    i_value = w1_value
                else:
                    i_value = w1_value / w2_value
                print >>fileD[op], "%d\t%f" % (i+1, i_value)
            elif op == 'average':
                i_value =(w1_value + w2_value) / 2
                print >>fileD[op], "%d\t%f" % (i+1, i_value)
            elif op == 'sum':
                i_value = w1_value + w2_value
                print >>fileD[op], "%d\t%f" % (i+1, i_value)
            #------------------------------------
        #--end op in opL--------------------------
    #------end i in posL-------------------------
#--------------------------------------------



def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    wig1 = options.wig1
    wig2 = options.wig2
    strand = options.strand
    global verbose
    verbose = options.verbose
    debug = options.debug
    outP = options.outP
    opL = options.op.split(',')
    fileD = {}
    for i in opL:
        fileD[i] = open(outP+'.'+i+'.wig', 'w')
    #-----------------------------------
    #opDict = {'d':mean, 'median':median, \
    #        'max':max, 'min':min, 'sum':sum}
    #--------------------------------
    outputChrList = set()
    outputSingle = 1
    wigDict1 = readWig(wig1, strand, outputSingle, outputChrList, 0,
            0)
    if verbose:
        print >>sys.stderr, outputChrList
    outputSingle = 0
    wigDict2 = readWig(wig2, strand, outputSingle, outputChrList,
            opL, fileD)
    outputChrList2 = list(outputChrList)
    lenocl = len(outputChrList2)
    wigDict = {}
    for i in range(lenocl):
        chr = outputChrList2[i]
        wigDict[chr] = {}
        wigDict[chr][0] = 0 # a random number
        processWig(wigDict, chr, outputChrList, opL, fileD)
    #print wigDict
    for i in fileD.values():
        i.close()
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



