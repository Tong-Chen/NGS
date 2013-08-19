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

File format: (for getShannon)
test_id 22    23    16    C   s
0610007C21Rik__XLOC_020042__=   7.14519 45.6536 46.4017 100.757 96.7285
0610007L01Rik__XLOC_020643__=   18.5009 7.2631  9.46834 11.1662 44.7087
0610007P08Rik__XLOC_006767__=   6.67136 5.67468 7.66136 4.43225 4.36321
0610007P14Rik__XLOC_006227__=   20.2475 26.611  40.441  58.4041 55.4284
0610007P22Rik__XLOC_011164__=   32.1198 55.6399 37.212  27.8336 18.7448
0610009B22Rik__XLOC_004475__=   8.67431 10.6717 13.0577 14.0596 22.4368
0610009D07Rik__XLOC_005329__=   46.6707 64.4448 70.8208 51.2847 43.1716
0610009O20Rik__XLOC_012558__=   8.14571 18.6023 17.8774 13.5515 13.2361


File format: (for getGene, normally the output from <getShannon>, 
The columns between the first and the last (not included) will be the
name of each lineage. They should be same ones as in sampleClass if
exists. Additional columns are unallowed and the last column is
shannon entropy.)
test_id	22	23	16	C	s	Shannon_index
0610007C21Rik__XLOC_020042__=	7.14519	45.6536	46.4017	100.757	96.7285	2.01988904654
0610007L01Rik__XLOC_020643__=	18.5009	7.2631	9.46834	11.1662	44.7087	1.97254835986
0610007P08Rik__XLOC_006767__=	6.67136	5.67468	7.66136	4.43225	4.36321	2.28663417741
0610007P14Rik__XLOC_006227__=	20.2475	26.611	40.441	58.4041	55.4284	2.21530259277
0610007P22Rik__XLOC_011164__=	32.1198	55.6399	37.212	27.8336	18.7448	2.23237916766
0610009B22Rik__XLOC_004475__=	8.67431	10.6717	13.0577	14.0596	22.4368	2.2429044139
0610009D07Rik__XLOC_005329__=	46.6707	64.4448	70.8208	51.2847	43.1716	2.29574228154
0610009O20Rik__XLOC_012558__=	8.14571	18.6023	17.8774	13.5515	13.2361	2.26868665144
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
from math import log as ln
from math import sqrt


def computeShannon(fh, header, verbose):
    for line in fh:
        if header:
            header -= 1
            print "%s\tShannon_index" % line.rstrip()
            continue
        #-------------------------
        lineL = line.split()
        expr = [float(i) for i in lineL[1:]]
        if verbose:
            print >>sys.stderr, expr
        expr_sum = sum(expr)
        if verbose:
            print >>sys.stderr, expr_sum
        assert expr_sum != 0
        expr_R = [1.0 * i / expr_sum for i in expr]
        if verbose:
            print >>sys.stderr, expr_R
        expr_Log = []
        for i in expr_R:
            if i != 0:
                expr_Log.append(i*ln(i)/ln(2))
            else:
                expr_Log.append(i)
        if verbose:
            print >>sys.stderr, expr_Log
        shannon = -1 * sum(expr_Log)
        print "%s\t%s" % (line.strip(),str(shannon))
    #-------------END reading file----------
#---------END computeShannnon----------------------

#def getGeneByShannon(fh,header,entropy,fpkm,cs):
def getSpecificByShannon(fh,header,entropy,fpkm,cs,strict_entropy):
    for line in fh:
        if header:
            header -= 1
            lineL = line.split()
            sampL = lineL[1:-1]
            if cs:
                for i in sampL:
                    assert i in cs.keys(), i
            lensampL = len(sampL)
            print "%s\tLineageInfo" % line.rstrip()
            continue
        #-------------------------
        lineL = line.split()
        expr = [float(i) for i in lineL[1:-1]]
        shannon = float(lineL[-1])
        lineageL = []
        if shannon <= entropy:
            for i92 in range(lensampL):
                if expr[i92] >= fpkm:
                    lineageL.append(sampL[i92])
            #---------------------------------------
            if 0 < len(lineageL) <= 2:
                print "%s\t%s" % (line.rstrip(),':'.join(lineageL))
            elif len(lineageL) > 2 and cs:
                tmp106 = set([cs[i] for i in lineageL])
                if len(tmp106) <= 2:
                    print "%s\t%s" % (line.rstrip(),':'.join(lineageL))
                #-----------strict_entropy--------------------------------
                elif strict_entropy != -1 and shannon < strict_entropy:
                    lineageL = []
                    mean, sd = mean_std(expr, lensampL)
                    #print >>sys.stderr, shannon, mean, sd, expr
                    large113 = mean + sd
                    for i114 in range(lensampL):
                        if expr[i114] > large113:
                            lineageL.append(sampL[i114])
                    #-----------------------------------------
                    if 0 < len(lineageL) <= 1:
                        print "%s\t%s" % (line.rstrip(),':'.join(lineageL))
                    elif len(lineageL) > 1 and cs:
                        tmp106 = set([cs[i] for i in lineageL])
                        if len(tmp106) <= 1:
                            print "%s\t%s" % (line.rstrip(),':'.join(lineageL))
            #-------------------------------------------
        #---------END each line------------
    #-------------END all lines---------------
#-----------END getGene-----------------------------------
def mean_std(data, len):
    sum_d = float(sum(data))
    mean_d = sum_d/len
    sd_d = sqrt(sum([(i-mean_d)**2 for i in data])*1.0/(len-1))
    return mean_d,sd_d
#-----------------------------------------------

def getCommonByShannon(fh,header,entropy,fpkm):
    for line in fh:
        if header:
            header -= 1
            lineL = line.split()
            sampL = lineL[1:-1]
            lensampL = len(sampL)
            print line,
            continue
        #-------------------------
        lineL = line.split()
        expr = [float(i) for i in lineL[1:-1]]
        shannon = float(lineL[-1])
        #lineageL = []
        lineage_cnt = 0
        if shannon > entropy:
            for i92 in range(lensampL):
                if expr[i92] >= fpkm:
                    #lineageL.append(sampL[i92])
                    lineage_cnt += 1
            #---------------------------------------
            if lineage_cnt == lensampL:
                print line,
        #-------One line------------------------
    #-----------All lines-----------------

#-------------------------------------------------
def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    desc = ""
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A tab separated file with the first \
column as names.")
    parser.add_option("-n", "--number_header", dest="header",
        metavar="a number", default=1, help="A number to indicate the header \
lines in the file. Default 1 means only one header line.")
    parser.add_option("-o", "--operation", dest="op",
        metavar="getShannon/getSpecific/getCommon", default='getShannon', help="Since this program is \
designed for two purposes, the operation you wanted should be chosen here. \
First <getShannon> : compute the shannon entropy for each gene. \
Second <getSpecific> : Get genes based on given criteria. \
Third <getCommon> : Get genes based on given criteria. \
Default <getShannon>.")
    parser.add_option("-e", "--entropy", dest="entropy",
        metavar="A number", default=2, help="The larger the entropy is, the \
more ubiquitous of a gene. Normally a small value would be better to \
get lineage specific genes and a large value would be better for \
housekeeping genes. It ranges from 0 to log2(N)[N is the number \
of lineages used here.] Default 2. If <getSpecific>, all genes with \
entropy no larger than than given value here will be processed furtherly. \
If <getCommon>, all genes with entropy larger than given value will \
processed furtherly.")
    parser.add_option("-E", "--strict-entropy", dest="strict_entropy",
        metavar="A number", default=-1, help="The paramter is only for \
<getSpecific>. This is used to extract some differential genes which \
expressed in all lineages but with much difference. Normally if an \
gene have FPKM larger than (mean+sd), it will be treated as specific. \
Default -1 means unsed. A number ranges from 0 to log2(N) is expected.")
    parser.add_option("-s", "--standard_for_high_expr", dest="fpkm",
        metavar="A number", default=1, help="This is mainly designed to \
filter genes by expression values.If a gene have rpkm \
larger than given value in one lineage and such high expression \
cannot be observed in more than two additional cell types,  it will \
be defined as a lineage specifc genes<for getSpecific>. If operation \
is <getCommon>, genes with entropy larger than given value here will \
be taken as common. Default 1.")
    parser.add_option("-m", "--mixed_sample", dest="sampleClass",
        metavar="", default='', help="Usually this parameter is \
useless unless you have very similar lineages. For example, in my \
sample, I have 'C1,C2,C3,E1,E2,F1'. I will assume \
'C1,C2,C3' may have similar expression patterns and 'E1,E2' \
may have similar expression patterns. Here I would give a tring like \
'C1:C2:C3;E1:E2;F1' to indicate my first three lineages can be \
treated as one lineage and the last two can be treated as one also. \
Make clear the usage of colon and semicolon and also the names given \
here must be the same as they appear in header line. The order \
dose not matter and single sample should be added too.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
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
    op = options.op
    entropy = float(options.entropy)
    strict_entropy = float(options.strict_entropy)
    header = int(options.header)
    verbose = options.verbose
    debug = options.debug
    fpkm = float(options.fpkm)
    cs = options.sampleClass
    if cs:
        tmpcs = [i.split(':') for i in cs.split(';')]
        a = len(tmpcs)
        csDict = {}
        for i in range(a):
            for j in tmpcs[i]:
                csDict[j] = i
        #--------------------------------------------------
        cs = csDict
        #print >>sys.stderr, cs
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    if op == 'getShannon':
        computeShannon(fh, header, verbose)
    elif op == 'getSpecific':
        getSpecificByShannon(fh,header,entropy,fpkm,cs,strict_entropy)   
    elif op == 'getCommon':
        getCommonByShannon(fh,header,entropy,fpkm)   
    else:
        print >>sys.stderr, "Unknow operation %s" % op
        sys.exit(1)
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



