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
    This is designed to do specified comparision from input file

Input file format:
    A   samp1
    B   samp1
    C   samp1
    B   samp2
    D   samp2
    A   samp3
    B   samp3

    OR
    A   samp1_up
    B   samp1_up
    C   samp1_up
    D   samp1_dw
    E   samp1_dw
    F   samp1_dw
    B   samp2_up
    D   samp2_dw
    A   samp3_up
    B   samp3_dw

Compare format:

    [samp1+samp2] --> Common items between samp1 and samp2
    [samp1+samp2+samp3] --> Common items among samp1, samp2 and samp3
    [samp1&samp2+samp3] --> Common items between combination of
                            (samp1, samp2) and samp3
    [samp1-samp2] --> samp1 special and samp2 special
    [samp1&samp2-samp3] --> union (samp1+samp2) special, and samp3
                            special
    
    When [up,up;dw,dw] is given:

    [samp1+samp2] --> Common items between 
                        samp1_up and samp2_up,
                        samp1_dw and samp1_dw
    [samp1+samp2+samp3] --> Common items among 
                        samp1_up, samp2_up and samp3_up
                        samp1_dw, samp2_dw and samp3_dw
    [samp1-samp2] --> Special items
                        samp1_up special and samp2_up special
                        samp1_dw special and samp2_dw special
    [samp1&samp2-samp3] --> Special items
                        union (samp1_up+samp2_up) special, samp3_up special
                        union (samp1_dw+samp2_dw) special, samp3_dw special

    When [up,dw;] is given (only for two sample comparison)

    [samp1+samp2] --> Common items between 
                        samp1_up and samp2_dw,
                        samp1_dw and samp1_up
    [samp1+samp2+samp3] --> Unsuitable
    [samp1-samp2] --> Special items
                        samp1_up special and samp2_dw special
                        samp1_dw special and samp2_up special
    [samp1&samp2-samp3] --> Special items
                        union (samp1_up+samp2_up) special, samp3_dw special
                        union (samp1_dw+samp2_dw) special, samp3_up special

'''
debug = 0

outputNameD = {}

import sys
import os
import re
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

def fprint(content):
    print >>sys.stderr, json_dumps(content,indent=1)

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
        metavar="FILEIN", help="A two column file with the first \
column containing gene names and the second column containing sample \
names. No header line is allowed.")
    parser.add_option("-o", "--out-prefix", dest="prefix",
        help="The output prefix of files.")
    parser.add_option("-s", "--sample-compare-patterns", dest="sample",
        help="The samples one want to compare. Usually in the format \
like (at most one <-> is allowed for each comparison, \
priority: & > + > - ; \
&: union; +: intersect; -: difference): \
'samp1+samp2;samp1+samp2+samp3;samp1-samp2;samp1&samp2-samp3;\
samp1&samp2-samp3&samp4'")
    parser.add_option("-S", "--sample-compare-name", dest="sampleName",
        help="Specify names for each sample comapre patterns. \
Optional. A list of ';' separated strings with the order as -s.")
    parser.add_option("-p", "--pair-mode", dest="pair_mode",
        help="Do complete pair comparsion. \
Accept '+' to get common items for each pair. \
Accept '-' to get special items for each pair. \
Accept '+-' for both.")
    parser.add_option("-P", "--pair-order", dest="pair_order",
        help="The order of comparison for each pair, like \
'a b c' will generate a-vs-b, a-vs-c, b-vs-c for '+' mode \
The order does not affect '-' mode, since special term for \
each will be output.")
    parser.add_option("-t", "--type", dest="type",
        help="Adding additional comparision by combining \
'.up,.up;.dw,.dw;.up,.dw' will generate \
<samp1.up, samp1.dw, samp2.up, samp2.dw> to specify additional layers. \
'.up' and '.up,.up' is the same. \
The dot (.) symbol can be changed based on the input file.")
    parser.add_option("-r", "--reverse", dest="reverse",
        default=1, help="Default 1 meaning for <-> two side \
special will be computed, accept 0 compute only one side special \
given in command line. Specially when <-S> is given, reverse \
will be assigned to 0.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def get_common(alist):
    '''
    alist = [set([1, 2, 3]), set([1, 2, 4]), ...]
    '''
    common = alist[0]
    for i in alist[1:]:
        common = common.intersection(i)
    return common

#--------get_common------------

def output(aDict, prefix):
    global outputNameD
    all = prefix + '.all'
    all = open(all, 'w')
    all_desp = prefix + '.all.desp'
    all_desp = open(all_desp, 'w')
    #print outputNameD
    #print aDict.keys()
    for label, valueS in aDict.items():
        label_sub = outputNameD.get(label, '')
        old_label = label
        if label_sub:
            label = label_sub
        file = prefix + '.' + label
        fh = open(file, 'w')
        print >>fh, '\n'.join(valueS)
        fh.close()
        for i in valueS:
            print >>all, "%s\t%s" % (i, label)
        print >>all_desp, "%s\t%s" % (old_label, label)
    #------------------
    all.close()
    all_desp.close()
#--------END output---------------
def GetIDs(labelL, dataD):
    '''
    labelL = ["samp1", ...] 
    '''
    if isinstance(labelL, str):
        return dataD.get(labelL, set())
    else:
        aSet = set()
        for label in labelL:
            aSet = aSet.union(dataD.get(label, set()))
        return aSet
#----------------------------
def intersect(compareSamp, op, dataD, aDict, reverse, sampName=''):
    '''
    compareSamp = [[samp1], [samp2]]
    compareSamp = [[samp1], [samp2], [samp3]]
    compareSamp = [[samp1, samp2], [samp3]]
    
    op = '+' or '-'

    dataD = [samp1: set([]), samp2: set([])]
    '''
    global outputNameD
    if debug:
        print >>sys.stderr, sampName
        print >>sys.stderr, outputNameD
    if debug:
        print >>sys.stderr, '#### intersect', compareSamp, op
    if op == '+':
        #newL = '__comm__'.join(['__union__'.join(j) \
        #    for i in compareSamp for j in i])
        newL = []
        for i in compareSamp:
            if isinstance(i, str):
                newL.append(i)
            elif isinstance(i, list):
                newL.append('__union__'.join(i))
        newL = '__comm__'.join(newL)
        newG = [GetIDs(i, dataD) for i in compareSamp]
        if debug:
            print >>sys.stderr, '##### newG', newG
            print >>sys.stderr, "##### dataD before common", dataD
        common = get_common(newG[:])
        aDict[newL] = list(common)
        outputNameD[newL] = sampName
        if debug:
            print >>sys.stderr, "##### dataD after common", dataD
            print >>sys.stderr, '###### Operation name', newL
            print >>sys.stderr, '###### Operation result', common
    elif op == '-':
        #firstL = '__comm__'.join(['__union__'.join(j) \
        #    for i in compareSamp[0] for j in i])
        firstL = []
        for i in compareSamp[0]:
            if isinstance(i, str):
                firstL.append(i)
            elif isinstance(i, list):
                firstL.append('__union__'.join(i))
        firstL = '__comm__'.join(firstL)
        secondL = []
        for i in compareSamp[1]:
            if isinstance(i, str):
                secondL.append(i)
            elif isinstance(i, list):
                secondL.append('__union__'.join(i))
        secondL = '__comm__'.join(secondL)

        #secondL = '__comm__'.join(['__union__'.join(j) \
        #    for i in compareSamp[1] for j in i])
        firstG = [GetIDs(i, dataD) for i in compareSamp[0]]
        secondG = [GetIDs(i, dataD) for i in compareSamp[1]]
        if debug:
            print >>sys.stderr, '##### firstG', firstG
            print >>sys.stderr, '##### secondG', secondG
            print >>sys.stderr, "##### dataD before common", dataD
        firstComm = get_common(firstG[:])
        secondComm = get_common(secondG[:])
        firstSpecial = firstComm.difference(secondComm)
        aDict['__special__'.join([firstL, secondL])] = list(firstSpecial)
        outputNameD['__special__'.join([firstL, secondL])] = sampName
        if reverse:
            secondSpecial = secondComm.difference(firstComm)
            aDict['__special__'.join([secondL, firstL])] = list(secondSpecial)
        if debug:
            print >>sys.stderr, "##### dataD after common", dataD
            print >>sys.stderr, '###### Operation name', \
            '__special__'.join([firstL,  secondL])
            print >>sys.stderr, '###### Operation result', firstSpecial
            if reverse:
                print >>sys.stderr, '###### Operation name', \
                '__special__'.join([secondL,  firstL])
                print >>sys.stderr, '###### Operation result', secondSpecial
    #---------------------------------------------
#-----------intersect--------------------------
def compare(compareSamp, typeL, op, dataD, aDict, reverse, sampName):
    if debug:
        print >>sys.stderr, "### compareSamp:", compareSamp, op
    for t in typeL:
        if len(t) == 2:
            assert len(compareSamp) == 2
            #compareSamp1 = [[j+t[0] for j in i] for i in compareSamp[0]]
            compareSamp1 = []
            for i in compareSamp[0]:
                if isinstance(i, str):
                    compareSamp1.append(i+t[0])
                elif isinstance(i, list):
                    compareSamp1.append([j+t[0] for j in i])
            #-------------------------------------
            compareSamp2 = []
            for i in compareSamp[1]:
                if isinstance(i, str):
                    compareSamp2.append(i+t[1])
                elif isinstance(i, list):
                    compareSamp2.append([j+t[1] for j in i])
            #-------------------------------------
            #compareSamp2 = [[j+t[1] for j in i] for i in compareSamp[1]]
            if debug:
                print >>sys.stderr, "### compareSamp1:", compareSamp1, op
                print >>sys.stderr, "### compareSamp2:", compareSamp2, op
            intersect([compareSamp1, compareSamp2], op, dataD, aDict,
                    reverse)
            compareSamp1 = [i+t[1] for i in compareSamp[0]]
            compareSamp2 = [i+t[0] for i in compareSamp[1]]
            if debug:
                print >>sys.stderr, "### compareSamp1:", compareSamp1, op
                print >>sys.stderr, "### compareSamp2:", compareSamp2, op
            intersect([compareSamp1, compareSamp2], op, dataD, aDict,
                    reverse)
        else:
            #compareSamp1 = [[[k+t[0] for k in j] for j in i] \
            #        for i in compareSamp]
            compareSamp1= []
            count = -1
            for i in compareSamp:
                count += 1
                compareSamp1.append([])
                for j in i:
                    if isinstance(j, str):
                        compareSamp1[count].append(j+t[0])
                    elif isinstance(j, list):
                        compareSamp1[count].append([k+t[0] for k in j])
                #-----------------------------
            if debug:
                print >>sys.stderr, "### compareSamp1:", compareSamp1, op
            intersect(compareSamp1, op, dataD, aDict, reverse,
                    sampName)
    #-----------------------Add additional later-
#----------Compare each---------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    prefix = options.prefix
    verbose = options.verbose
    global debug
    debug = options.debug
    reverse = int(options.reverse)
    #------------Get sampD----------------------
    # samp1+samp2;samp1+samp2+samp3;
    # samp1&samp2+samp3; different when up;dw is used
    # samp1-samp2;
    # samp1&samp2-samp3;samp1&samp2-samp3&samp4'
    samp = options.sample
    sampName = options.sampleName
    sampNameD = {'+':[], '-':[]}
    if sampName:
        sampNameL = sampName.split(';')
        reverse = 0
        if debug:
            print >>sys.stderr, sampNameL
    '''
    sampD = {'+': [
                [[samp1], [samp2]],
                [[samp1], [samp2], [samp3]],
                [[samp1, samp2], samp3],
                ], 
             '-': [
                [[samp1], [samp2]], 
                [[samp1, samp2], [samp3]], 
                [[samp1, samp2], [samp3, samp4]], 
                ]}
    '''
    sampD = {'+': [], '-': []}
    if samp:
        sampL = samp.split(';')
        i348 = -1
        for samp in sampL:
            i348 += 1
            combine  = samp.count('+')
            subtract = samp.count('-')
            #assert combine * subtract == 0, "Illegel samp: %s" % samp
            assert subtract <= 1, "Illegel samp: %s" % samp
            #if subtract:
            #    for twoDiff in samp.split('-'):

            if subtract:
                sampD['-'].append([[j.split('&') for j in i.split('+')] \
                    for i in samp.split('-')])
                if sampName:
                    sampNameD['-'].append(sampNameL[i348])
            elif combine:
                sampD['+'].append([i.split('&') for i in samp.split('+')])
                if sampName:
                    sampNameD['+'].append(sampNameL[i348])
            #--------------------------------------
        #-----------------------------------
        if debug:
            fprint(sampD)
            print >>sys.stderr, sampNameD
    #------------Get sampD----------------------
    #--Get pair---------------------------------
    samp_pair = options.pair_mode
    pair_order = options.pair_order

    type = options.type
    '''
    type = [['up'], ['dw'], ['up', 'dw']]
    '''
    if type:
        typeL = [list(set(i.split(','))) for i in type.split(';')] 
    else:
        typeL = [['']]
    if debug:
        print >>sys.stderr, "## typeL:", typeL
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    dataD = {}
    aDict = {}
    for line in fh:
        value, key = line.strip().split('\t')
        if key not in dataD:
            dataD[key] = set()
        dataD[key].add(value)
    #-------------END reading file----------
    if debug:
        print "# dataD", dataD
    #----Combine-----------------------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------

    '''
    sampD = {'+': [
                [[samp1], [samp2]],
                [[samp1], [samp2], [samp3]],
                [[samp1, samp2], samp3],
                ], 
             '-': [
                [[samp1], [samp2]], 
                [[samp1, samp2], [samp3]], 
                [[samp1, samp2], [samp3, samp4]], 
                ]}
    '''
    for op in '+-': 
        if op in sampD:
            sampL = sampD[op]
            sampNameL = sampNameD.get(op)
            i425 = -1
            for compareSamp in sampL:
                i425 += 1
                if sampNameL:
                    sampName = sampNameL[i425]
                else:
                    sampName = ''
                compare(compareSamp, typeL, op, dataD, aDict, reverse,
                        sampName)
        #-------Compare all sampD----------------
    #--------Compare all sampD-----------------------
    if samp_pair:
        if pair_order:
            pair_orderL = pair_order.split(' ')
        else:
            pair_orderL = dataD.keys()
            pair_orderL.sort()
        len1 = len(pair_orderL)
        for i in range(len1-1):
            samp1 = [pair_orderL[i]]
            for j in range(i+1, len1):
                samp2 = [pair_orderL[j]]
                compareSamp = [samp1, samp2]
                if debug:
                    print "## compareSamp", compareSamp
                for op in samp_pair:
                    compare(compareSamp, typeL, op, dataD, aDict)
        #----Generate compare-pair---
    #----------END samp_pair---------------
               
    if debug:
        fprint(aDict)
    output(aDict, prefix)
                    

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


#    for op in '+-': 
#        if op in sampD:
#            sampL = sampD[op]
#            for compareSamp in sampL:
#                compare(compareSamp, typeL, op, dataD, aDict)
#                for t in typeL:
#                    if len(t) == 2:
#                        assert len(compareSamp) == 2
#                        compareSamp1 = [i+'.'+t[0] for i in compareSamp[0]]
#                        compareSamp2 = [i+'.'+t[1] for i in compareSamp[1]]
#                        aDict = intersect([compareSamp1, compareSamp2], op, dataD)
#                        output(aDict, prefix)
#                        compareSamp1 = [i+'.'+t[1] for i in compareSamp[0]]
#                        compareSamp2 = [i+'.'+t[0] for i in compareSamp[1]]
#                        aDict = intersect([compareSamp1, compareSamp2], op, dataD)
#                        output(aDict, prefix)
#                    else:
#                        compareSamp1 = [[j+'.'+t[0] for j in i] for i in compareSamp]
#                        aDict = intersect(compareSamp1, op, dataD)
#                        output(aDict, prefix)
#                #-----------------------Add additional later-
#            #----------Compare each---------------------
#        #-------Compare all sampD----------------
#    #--------Compare all sampD-----------------------
