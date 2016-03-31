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
desc = '''
Functional description:
    This is designed to summary the enrichment of peaks along
    speicified regions.

Input file:

1. The simplest

344@given.1 1
344@given.2 1
344@up.1 1
344@up.2 1
344@dw.1 1
344@dw.2 1

2. File with coordination information

chr1	1	20	344@given.1 1
chr1	1	20	344@given.2 1
chr1	1	20	344@up.1 1
chr1	1	20	344@up.2 1
chr1	1	20	344@dw.1 1
chr1	1	20	344@dw.2 1

3. File with complex score (normally this is the output of groupBy)

chr2 54730000 54730500 344@given.1 CTCF:1,FOXM1:1,RAD21:1,RUNX3:1,SMC3:
chr2 54730500 54731000 344@given.2 CTCF:1,FOXM1:1,RAD21:1,RUNX3:1,SMC3:
chr2 54731000 54731500 344@given.3 .:1,
chr2 54731500 54732000 344@given.4 .:1,

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

debug = 0

def fprint(content):
    print json_dumps(content,indent=1)

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="Input file formated as exmples.")
    parser.add_option("-l", "--label", dest="labelC",
        default=4, help="The column containing the names of each \
position (1-based). Default 4.")
    parser.add_option("-L", "--labelType", dest="labelType",
        default="@-2", help="The default parameter is <@-2> indicating \
the real label can be got by separating the label column using <@> \
then extracing the second element (1-based).")
    parser.add_option("-o", "--labelOrder", dest="labelOrder",
        default="up.40.given.20.dw.40", 
        help="The string can be used to generate the full list of \
labels. Default <up.40.given.20.dw.40> will generate \
<up.1 up.2 ... up.40 given.1 given.2 ... given.20 dw.1 dw.2 ... dw.40>")
    parser.add_option("-T", "--labelNumType", dest="labelNumType",
        default=0, 
        help="Transfer labels to numbers. Eg, \
When <labelOrder> is <up.40.given.20.dw.40>, true labels would be \
<up.1 up.2 ... up.40 given.1 given.2 ... given.20 dw.1 dw.2 ... dw.40>. \
When <labelNumType> is <1>, true labels wil be transfered into \
<1 2 ... 40 41 42 ... 60 61 62 ... 100>. Default <0> indicating \
no transfer,  accept <1> to transfer.")
    parser.add_option("-s", "--score", dest="scoreC",
        default=5, help="The column containing the scores of each \
position (1-based). Default 5.")
    parser.add_option("-t", "--score-type", dest="scoreT",
        default='complex', help="Specify the type of score column. \
Currently two types are supported, when there is only one number in \
score column, <simple> should be given here. <complex> is the default \
parameter to deal with the score column indicated in example data.")
    parser.add_option("-N", "--rename-score", dest="scoreN",
        default='score', help="When <--score-type> is <simple>, \
a string given here will be used to represent the name of score \
column. Default <score>. Normally sample information should be \
given here.")
    parser.add_option("-n", "--normalize", dest="normalize",
        default='2-3', help="Specify one column containing \
normalization factor (the length of the region) or two columns \
containing the <start-end>(default 2-3) columns for computing region length. \
<FALSE> indicating no normalize.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def normalizeScore(score, lineL, normalize):
    if normalize == 'FALSE':
        return score
    else:
        normalizeL = [int(i)-1 for i in normalize.split('-')]   
        lennorm = len(normalizeL)
        if lennorm == 1:
            region = float(lineL[normalizeL[0]])
        elif lennorm == 2:
            region = float(lineL[normalizeL[1]]) - float(lineL[normalizeL[0]])
        return score / region
#----------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    labelC = int(options.labelC) - 1
    #@-2
    labelType = options.labelType
    labelSep, labelGet = labelType.split('-')
    labelGet = int(labelGet)-1
    #up.50.given.20.dw.50
    labelOrder = options.labelOrder
    labelOrderL = labelOrder.split('.')
    labelL = []
    summary = 1
    for i in range(1,len(labelOrderL),2):
        count = int(labelOrderL[i])
        summary += count
        type = labelOrderL[i-1]
        tmpL = [type+'.'+str(j+1) for j in range(count)]
        labelL.extend(tmpL)
    
    labelNumType = int(options.labelNumType)
    labelNumTypeD = {}
    if labelNumType:
        for i in range(1, summary):
            labelNumTypeD[labelL[i-1]] = str(i)
    #-----------------------------------------
    scoreC = int(options.scoreC) - 1
    scoreN = options.scoreN
    #simple complex
    scoreType = options.scoreT
    #FALSE 2 2-3
    normalize = options.normalize
    #--------------------------------------------------------
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    '''
    aDict = {'given.1':{'CTCF':2, 'RAD21':3}, 
            'given.2':{'CTCF':2, 'RAD21':3}, 
            }

    aDict = {'given.1':{scoreN:2}, 
            'given.2':{scoreN:2}, 
            }
    '''
    aDict = {}
    aSet = set()
    for line in fh:
        lineL = line.strip().split('\t')
        label = lineL[labelC].split(labelSep)[labelGet]
        score = lineL[scoreC]
        if label not in aDict:
            aDict[label] = {}
        if scoreType == 'simple':
            if scoreN not in aDict[label]:
                aDict[label][scoreN] = []
            score = normalizeScore(float(score), lineL, normalize)
            aDict[label][scoreN].append(score)
        elif scoreType == 'complex':
            #print score
            scoreL = [i.split(':') for i in score.rstrip(',').split(',')]
            #print scoreL
            for scoreEle in scoreL:
                peak_name, peak_value = scoreEle
                if peak_name == '.':
                    continue
                peak_value = normalizeScore(float(peak_value), lineL, normalize)
                if peak_name not in aDict[label]:
                    aDict[label][peak_name] = []
                aDict[label][peak_name].append(peak_value)
        #-----------------------------------------
    #-------------END reading file----------
    bDict = {}
    '''
    bDict = {'CTCF':{'given.1':1, 'given.2':2},  
            'RAD21':{'given.1':1, 'given.2':2}}

    bDict = {scoreN:{'given.1':1, 'given.2':2}}
    '''
    for label, labelD in aDict.items():
        for peak, peak_valueL in labelD.items():
            peak_value = sum(peak_valueL)
            if peak not in bDict:
                bDict[peak] = {}
            bDict[peak][label] = str(peak_value)
    #-----------output-------------------------
    print 'Pos\tType\tValue'
    for peak, labelD in bDict.items():
        print '\n'.join(['\t'.join([labelNumTypeD.get(label, label), 
            peak, labelD.get(label,'0')]) \
            for label in labelL])
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
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
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


