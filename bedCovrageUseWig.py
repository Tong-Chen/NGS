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
'''

import collections
from numpy import mean,median,max,min,sum
from numpy import array as np_array
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
        metavar="BED_REGION", help="Regions in bed file format, SORT \
BY chromsome.***")
    parser.add_option("-w", "--wig", dest="wig",
        metavar="WIG", help="Regions in wig file format.***")
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

def readWig(wig, strand):
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
            pos_fixed = 0
            chromi = i.rfind('chrom=')
            assert chromi != -1, "Wrong format no chr %s" % i
            chr = i[chromi+6:].strip().split()[0]
            spani = i.rfind("span=")
            span = 1 if spani == -1 else \
                int(i[spani+5:].strip().split()[0])
            pos = 1
            neg = 2
        elif i.startswith("fixedStep"):
            chromi = i.rfind('chrom=')
            assert chromi != -1, "Wrong format no chr %s" % i
            chr = i[chromi+6:].strip().split()[0]
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
            for position in xrange(start, start+span):
                if not strand:
                    wigDict[chr][position] = float(lineL[pos])
                else:
                    wigDict[chr][position] = array('f', \
                        [float(lineL[pos]), float(lineL[neg])])
                #--------------------------------------
            #-----------------------------------------
        #--------end processing one line ----------------------
    #----------------END reading whole file-------------------
    return wigDict
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


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    bed = options.bed
    wig = options.wig
    strand = options.strand
    verbose = options.verbose
    debug = options.debug
    wigDict = readWig(wig, strand)
    opL = options.op.split(',')
    name_mode = int(options.name)
    #-----------------------------------
    opDict = {'mean':mean, 'median':median, \
            'max':max, 'min':min, 'sum':sum}
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(bed)
    #--------------------------------
    if name_mode:
        print "#name\t%s" % '\t'.join(opL)
    else:
        print "#%s" % '\t'.join(opL)
    for line in fh:
        lineL = line.strip().split('\t')
        chr   = lineL[0]
        start = int(lineL[1])
        end   = int(lineL[2])
        innerD = wigDict[chr]
        if strand:
            strand_in = lineL[5]
            if name_mode:
                name = lineL[3]+'@'+ strand_in
            else:
                name = line.strip()
            strand_num = 0 if strand_in=='+' else 1
            valueL = np_array([innerD.get(i,[0,0])[strand_num] \
                for i in xrange(start,end)])
        else:
            if name_mode:
                name = lineL[3]
            else:
                name = line.strip()
            valueL = np_array([innerD.get(i,0) \
                for i in xrange(start,end)])
        #----------------------------------------------
        #print valueL
        #for op in opL:
        #    tmpL.append(str(opDict[op](valueL)))
        print "%s\t%s" % (name, \
            '\t'.join([str(opDict[op](valueL)) for op in opL]))
    #-------------END reading file----------
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



