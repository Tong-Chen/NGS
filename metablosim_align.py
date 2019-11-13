#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to align metabolism data.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
from bisect import bisect_left
from fromMassToStructure import binarySearch, findTheNearest
#from multiprocessing.dummy import Pool as ThreadPool
from airflowDagGen import iterateParameters2

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

def fprint(content):
    """ 
    This is a Google style docs.

    Args:
        param1(str): this is the first param
        param2(int, optional): this is a second param
            
    Returns:
        bool: This is a description of what is returned
            
    Raises:
        KeyError: raises an exception))
    """
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
    parser.add_option("-i", "--main-file", dest="filein",
        metavar="FILEIN", help="Main metabolism file.")
    parser.add_option("-a", "--other-files", dest="other",
        metavar="FILEIN", help="A list of files waiting to be aligned to main file.")
    parser.add_option("-l", "--other-files-labels", dest="other_labels",
        help="A list of strings withs same order as those given to <-a> to label each file.")
    parser.add_option("-b", "--base-peak-retention-time", dest="base_peak_rt",
        help="Base peak retention time for each file separated by ';'")
    parser.add_option("-c", "--common-peak-retention-time", dest="common_peak_rt",
        help="Common peak retention time for each file separated by ';'")
    parser.add_option("-t", "--retention-time-sd", 
        dest="retention_time_sd", 
        default=0.01, type='float', 
        help="Allowed deviation of retention time. Default 0.01.")
    parser.add_option("-m", "--mz-ratio-sd", 
        dest="mz_ratio_sd", default=0.01, type='float', 
        help="Allowed deviation of M/Z. Default 0.01.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readInfile(file):
    header = 1
    '''
    Compound_name   RT  m/z
    10.19_227.1274m/z   10.19   227.1274385
    11.06_717.2171m/z   11.06   717.2170699
    6.47_303.1798m/z    6.47    303.1797728
    3.36_569.1821m/z    3.36    569.1820643
    '''
    mzD = {}
    '''
    mzD = {mz1:[[label1, rt1], [label2, rt2]], mz2:[[label3, rt3]]}
    '''
    for line in open(file):
        if header:
            header -= 1
            continue
        lineL = line.strip().split('\t')
        label = lineL[1]+'_'+lineL[2]
        #label, rt, mz, null = line.split('\t', 3)
        rt = float(lineL[1])
        mz = float(lineL[2])
        if mz not in mzD:
            mzD[mz] = []
        mzD[mz].append([label, rt])
    return mzD
#---------readInfile------------------------

def alignTomain(primaryMZd, primaryMZd_keys, secondaryMZd, primary_base_rt, primary_common_rt, sec_base_rt, sec_common_rt, retention_time_sd, mz_ratio_sd, matchD, matchE, label):
    '''
    primaryMZd = = {mz1:[[label1, rt1], [label2, rt2]], mz2:[[label3, rt3]]}
    matchD = {'primarymz_label': {samp_label1:[
                                    [semzLabel1, secmz1, sec_rt1], 
                                    [semzLabel2, secmz2, sec_rt2], 
                                    ..], 
                                  samp_label2:[forth_mzL]}, 
              'primarymz_label2':{samp_label3:[forth_mzL]}}
    '''
    secondaryMZd_keys = secondaryMZd.keys()
    secondaryMZd_keys.sort()
    sec2priD = binarySearch(primaryMZd_keys, secondaryMZd_keys, mz_ratio_sd)
    '''
    sec2priD = {secmz1: [primarymz1, primatymz2, ...], secmz2: [], secmz3:[primarymz3]}
    '''
    pri_c_b_rt = primary_common_rt-primary_base_rt
    sec_c_b_rt = sec_common_rt-sec_base_rt
    if debug:
        print >>sys.stderr, "Input parameters:"
        print >>sys.stderr, "primary_common_rt={}; primary_base_rt={};".format(primary_common_rt, primary_base_rt)
        print >>sys.stderr, "sec_common_rt={}; sec_base_rt={};".format(sec_common_rt, sec_base_rt)
        print >>sys.stderr, "pri_c_b_rt={}; sec_c_b_rt={}".format(pri_c_b_rt, sec_c_b_rt)
        print >>sys.stderr, '-------------------------'
    for secmz, primarymzL in sec2priD.items():
        #if secmz not in matchD:
        #    matchD[secmz] = {}
        if debug:
            print >>sys.stderr, '\n-------------------------'
            print >>sys.stderr, "+Search MZ: {}".format(secmz)
            print >>sys.stderr, 'Target MZ:'+ ' ' + ', '.join([str(i) for i in primarymzL])
        for sec_label, sec_rt in secondaryMZd[secmz]:
            #matchD[secmz][sec_label] = sec_label
            sec_rt_ratio = (sec_rt-sec_base_rt) / sec_c_b_rt
            if debug:
                print >>sys.stderr, '\nsec_rt={};sec_rt_ratio={}'.format(sec_rt, sec_rt_ratio)
            matched = 0
            transferedRT = pri_c_b_rt/sec_c_b_rt*(sec_rt-sec_base_rt)+primary_base_rt
            if primarymzL:
                #for primarymz in primarymzL:
                #    print >>sys.stderr, primaryMZd[primarymz]
                primary_rt_ratioL = [(eachp[1]-primary_base_rt)/pri_c_b_rt for primarymz in primarymzL for eachp in primaryMZd[primarymz]]
                pos_tmp = 0
                # Add position to each value for index
                primary_rt_ratioL_order = []
                for i in primary_rt_ratioL:
                    primary_rt_ratioL_order.append([i, pos_tmp])
                    pos_tmp += 1
                primary_rt_ratioL_order.sort(key=lambda x: x[0])
                primary_rt_ratioL = [i[0] for i in primary_rt_ratioL_order]
                if debug:
                    print >>sys.stderr, 'RT ratio: '+','.join([str(i) for i in primary_rt_ratioL])
                    print >>sys.stderr, primary_rt_ratioL_order
                    print >>sys.stderr, 'Original RT: '+','.join([str(eachp[1]) for primarymz in primarymzL for eachp in primaryMZd[primarymz]])
                primary_label_tmpL = [eachp[0] for primarymz in primarymzL for eachp in primaryMZd[primarymz]]
                primary_rt_ratioL.sort()
                closest_value, closest_pos = findTheNearest(sec_rt_ratio, primary_rt_ratioL, 0, len(primary_rt_ratioL))
                closest_pos = primary_rt_ratioL_order[closest_pos][1]
                if debug:
                    print >>sys.stderr, "closest_value={}; closest_pos={}".format(closest_value, closest_pos)
                retention_time_compute = abs(closest_value - sec_rt_ratio)
                if retention_time_compute <= retention_time_sd:
                    primarymz_label = primary_label_tmpL[closest_pos]   
                    if primarymz_label not in matchD:
                        matchD[primarymz_label] = {}
                    if label not in matchD[primarymz_label]:
                        matchD[primarymz_label][label] = []
                    # sec_label: label for secondary metabolism
                    # 
                    matchD[primarymz_label][label].append([sec_label,secmz, retention_time_compute, transferedRT])
                    #matchD[primarymz_label][label].append(sec_label)
                    matched = 1
                    if debug:
                        print >>sys.stderr, "@Find match: {}-->{}".format(primarymz_label, sec_label)
            if not matched:
                transferedRT = pri_c_b_rt/sec_c_b_rt*(sec_rt-sec_base_rt)+primary_base_rt
                if secmz not in primaryMZd:
                    primaryMZd[secmz] = []
                new_primary_label = str(transferedRT)+'_'+str(secmz)
                primaryMZd[secmz].append([new_primary_label, transferedRT])
                if new_primary_label not in matchD:
                    matchD[new_primary_label] = {}
                if label not in matchD[new_primary_label]:
                    matchD[new_primary_label][label] = []
                matchD[new_primary_label][label].append([sec_label])
                if debug:
                    print >>sys.stderr, "@No match: {}-->{}".format(new_primary_label, sec_label)
            #----------------------------
        #-------------------Iterate each label have same mz in secondary file
    #--------------------------------
    """
    matchD = {'primarymz_label': {samp_label1:[
                                  [semzLabel1,secmz1,retention_time_compute,transferedRT], 
                                  [semzLabel2,secmz2,retention_time_compute,transferedRT], 
                                    ..], 
    """
    #matchE = {}
    for primarymz_label, samplabelD in matchD.items():
        if primarymz_label not in matchE:
            matchE[primarymz_label] = {}
        else:
            print >>sys.stderr, "Duplicate *"+primarymz_label+'*'
            sys.exit(1)
        #-------------------------
        for label, match_209L in samplabelD.items():
            if label not in matchE[primarymz_label]:
                matchE[primarymz_label][label] = []
            else:
                print >>sys.stderr, "Duplicate *"+label+'*'
                sys.exit(1)
            #----------------------------------------------------
            if (len(match_209L) > 1):
                print >>sys.stderr, "****"
                print >>sys.stderr, match_209L
                match_209L.sort(key=lambda x: x[2])
                matchE[primarymz_label][label].append(match_209L[0][0])
                print >>sys.stderr, "***"
                print >>sys.stderr, match_209L[1:]
                for match_209 in match_209L[1:]:
                    print >>sys.stderr, "**"
                    print >>sys.stderr, match_209
                    sec_label, secmz, retention_time_compute, transferedRT = match_209 
                    if secmz not in primaryMZd:
                        primaryMZd[secmz] = []
                    new_primary_label = str(transferedRT)+'_'+str(secmz)
                    primaryMZd[secmz].append([new_primary_label, transferedRT])
                    if new_primary_label not in matchE:
                        matchE[new_primary_label] = {}
                    if label not in matchE[new_primary_label]:
                        matchE[new_primary_label][label] = []
                    matchE[new_primary_label][label].append(sec_label)
                    #if debug:
                    #    print >>sys.stderr, "@No match: {}-->{}".format(new_primary_label, sec_label)
            else:
                matchE[primarymz_label][label].append(match_209L[0][0])
                #pass
       # for label in secondaryLabelL:
       #     match_labelL = matchD_pri_sec.get(label, ['NA'])
       #     match_labelL.sort(key=lambda x: x[2])

    #return matchE
#-----END alignTomain-----------------



#-------------------------alignTomain----------------
def outputMatchship(primaryMZd_labels, matchD, secondaryLabelL):
    '''
    matchD = {'primarymz_label': {samp_label1:[
                                    semzLabel1, 
                                    semzLabel2, 
                                    ..], 
    #matchD = {'primarymz_label': {samp_label1:[secondary_mzL, third_mzL, ..], 
    #                              samp_label2:[forth_mzL]}, 
    #          'primarymz_label2':{samp_label3:[forth_mzL]}}
    '''
    print '{}\t{}'.format("Primary", '\t'.join(secondaryLabelL))
    count = 0
    for i in primaryMZd_labels:
        if 0:
            print >>sys.stderr, i
            print >>sys.stderr, matchD.get(i)
            count += 1
            if count > 100:
                sys.exit(1)
        else:
            tmpL = [';'.join(matchD.get(i, {}).get(label, ['NA'])) for label in secondaryLabelL]
            print '{}\t{}'.format(i, '\t'.join(tmpL))
#---------------------------------------

def outputMatchship_unfinished(primaryMZd_labels, matchD, secondaryLabelL):
    '''
    matchD = {'primarymz_label': {samp_label1:[
                                    [semzLabel1, secmz1, sec_rt1], 
                                    [semzLabel2, secmz2, sec_rt2], 
                                    ..], 
    #matchD = {'primarymz_label': {samp_label1:[secondary_mzL, third_mzL, ..], 
    #                              samp_label2:[forth_mzL]}, 
    #          'primarymz_label2':{samp_label3:[forth_mzL]}}
    '''
    print '{}\t{}'.format("Primary", '\t'.join(secondaryLabelL))
    for i in primaryMZd_labels:
        matchD_pri_sec = matchD.get(i, {})
        for label in secondaryLabelL:
            match_labelL = matchD_pri_sec.get(label, ['NA'])
            if len(match_labelL) == 1:
                print '{}\t{}'.format(i, '\t'.join(tmpL))
            else:
                match_labelL.sort(key=lambda x: x[2])


        #tmpL = [';'.join(matchD.get(i, {}).get(label, ['NA'])) for label in secondaryLabelL]
        #print '{}\t{}'.format(i, '\t'.join(tmpL))
#---------------------------------------


def outputMatchship_multiple(primaryMZd_labels, matchD, secondaryLabelL):
    '''
    matchD = {'primarymz_label': {samp_label1:[secondary_mzL, third_mzL, ..], 
                                  samp_label2:[forth_mzL]}, 
              'primarymz_label2':{samp_label3:[forth_mzL]}}
    '''
    print '{}\t{}'.format("Primary", '\t'.join(secondaryLabelL))
    for i in primaryMZd_labels:
        tmpL = [';'.join(matchD.get(i, {}).get(label, ['NA'])) for label in secondaryLabelL]
        #tmpL2 = [i.split(';') for i in tmpL]
        tmpL = '{}\t{}'.format(i, '\t'.join(tmpL))
        tmpL = [j.split(';') for j in tmpL.split('\t')]
        resultL = []
        recordL = []
        iterateParameters2(tmpL, recordL, resultL)
        print '\n'.join(['\t'.join(i) for i in resultL])
#---------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    secondaryL = options.other.split()
    secondaryLabelL = options.other_labels.split()
    base_peak_rtL = [float(i) for i in options.base_peak_rt.split(';')]
    common_peak_rtL = [float(i) for i in options.common_peak_rt.split(';')]
    retention_time_sd = options.retention_time_sd
    mz_ratio_sd = options.mz_ratio_sd
    
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    matchD = {}
    matchE = {}
    '''
    matchD = {'primarymz_label': []}
    '''
    primaryMZd = readInfile(file)
    #primaryMZd_keys = primaryMZd.keys()
    #primaryMZd_keys.sort()
    #primary_base_peak = base_peakL[0]
    #primary_common_peak = common_peakL[0]
    primary_base_rt = base_peak_rtL[0]
    primary_common_rt = common_peak_rtL[0]
    #secondaryMZd = {}
    for label, secondary,sec_base_rt, sec_common_rt in zip(secondaryLabelL, secondaryL, base_peak_rtL[1:], common_peak_rtL[1:]):
        primaryMZd_keys = primaryMZd.keys()
        primaryMZd_keys.sort()
        secondaryMZd = readInfile(secondary)
        #sec_base_rt = secondaryMZd[sec_base_peak]
        #sec_common_rt = secondaryMZd[sec_common_peak]
        alignTomain(primaryMZd, primaryMZd_keys, secondaryMZd, primary_base_rt, primary_common_rt, sec_base_rt, sec_common_rt, retention_time_sd, mz_ratio_sd, matchD, matchE, label)
        
    primaryMZd_keys = primaryMZd.keys()
    '''
    mzD = {mz1:[[label1, rt1], [label2, rt2]], mz2:[[label3, rt3]]}
    '''
    primaryMZd_labels = [j[0] for i in primaryMZd.values() for j in i]

    """"
    matchE = {'primarymz_label': {samp_label1:[semzLabel1], 
    """
    outputMatchship(primaryMZd_labels, matchE,secondaryLabelL)
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


