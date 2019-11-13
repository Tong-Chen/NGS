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

    1. Transpose data to make each row is one sample.
    2. Uppercase sample names
    3. Match sample names to normal name (tcga-5l-aat0 to TCGA-5L-AAT0-01)
    4. Correct survival data for survival analysis (dead:1; unknown: 0) [vital_status, days_to_death, days_to_last_followup, days_to_last_known_alive will be ued]

Clinical data file:

Hybridization REF	tcga-5l-aat0	tcga-5l-aat1	tcga-a1-a0sp	tcga-a2-a04v
Composite Element REF	value	value	value	value
years_to_birth	42	63	40	39
vital_status	0	0	0	1
days_to_death	NA	NA	NA	1920
days_to_last_followup	1477	1471	584	NA
tumor_tissue_site	breast	breast	breast	breast
pathologic_stage	stage iia	stage iv	stage iia	stage iia
pathology_T_stage	t2	t2	t2	t2
pathology_N_stage	n0	n0	n0 (i-)	n0 (i-)
pathology_M_stage	m0	m1	m0	m0
gender	female	female	female	female
date_of_initial_pathologic_diagnosis	2010	2010	2007	2005
days_to_last_known_alive	NA	NA	NA	NA
radiation_therapy	yes	no	NA	no
histological_type	infiltrating lobular carcinoma	infiltrating lobular carcinoma	infiltrating ductal carcinoma	infiltrating ductal carcinoma
number_of_lymph_nodes	0	0	0	0
race	white	white	NA	white
ethnicity	hispanic or latino	hispanic or latino	not hispanic or latino	not hispanic or latino

sampleID file

TCGA: Project
A1  : Tissue source site
A0SB: Participant
01-09: Tumor types
10-19: normal types
20-29: control types

HYBRIDIZATION R
TCGA-A1-A0SB-01
TCGA-A1-A0SD-01
TCGA-A1-A0SE-01
TCGA-A1-A0SF-01
TCGA-A1-A0SG-01
TCGA-A1-A0SH-01
TCGA-A1-A0SI-01
TCGA-A1-A0SJ-01
TCGA-A1-A0SK-01

Attribute file

HYBRIDIZATION R	ATTR1
TCGA-A1-A0SB-01	1
TCGA-A1-A0SD-01	1
TCGA-A1-A0SE-01	1
TCGA-A1-A0SF-01	1
TCGA-A1-A0SG-01	1
TCGA-A1-A0SH-01	1
TCGA-A1-A0SI-01	1
TCGA-A1-A0SJ-01	1
TCGA-A1-A0SK-01	1


'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0


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
        metavar="FILEIN", help="Input firehose cilincal data")
    parser.add_option("-s", "--sample-id-file", dest="sample_id",
        help="Optional. With content listed above. The program will try to match original id with this updated IDs if given.")
    parser.add_option("-a", "--additional-file", dest="additional_attr",
        help="Optional. With first column as sample id and other columns as sample attributes.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def getSurvialData(aDict):
    """
    vital_status, days_to_death, days_to_last_followup, days_to_last_known_alive will be used

    Days.survial column will be added.
    """
    for sample, sampleD in aDict.items():
        assert 'Days.survial' not in sampleD
        survial = sampleD['days_to_death']
        if survial == 'NA':
            survial = sampleD['days_to_last_followup']
        if survial == 'NA':
            survial = sampleD['days_to_last_known_alive']
        #assert survial != 'NA', 'must have survial value not NA for '+sample
        sampleD['Days.survial'] = survial

#-------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    sample_id = options.sample_id
    verbose = options.verbose
    additional_attr = options.additional_attr
    global debug
    debug = options.debug
    #-----------------------------------
    sampleIDd = {}
    typeD = {}
    if sample_id:
        for line in open(sample_id):
            if line.find('-') != -1:
                sample_full = line.strip()
                sample_fullL = sample_full.split('-')
                sample_pre = '-'.join(sample_fullL[:3])
                type = sample_fullL[3]
                #01-09: Tumor types
                #10-19: normal types
                #20-29: control types
                if type[0] == '0':
                    typeD[sample_full] = 'Tumor_type'
                elif type[0] == '1':
                    typeD[sample_full] = 'Normal_type'
                elif type[0] == '2':
                    typeD[sample_full] = 'Control_sample'

                #assert sample_pre not in sampleIDd, "Duplicate "+sample_full+'. STH WRONG!!!'
                if sample_pre not in sampleIDd:
                    sampleIDd[sample_pre] = [sample_full]
                else:
                    sampleIDd[sample_pre].append(sample_full)
                #------------------------------------------
        if debug:
            print >>sys.stderr, sampleIDd
    #-------------------------------------------------------------
    attrD = {}
    if additional_attr:
        header = 1
        for line in open(additional_attr):
            lineL = line.strip().split('\t')
            if header:
                attributeL = lineL[1:]
                header -= 1
                continue
            #------------------------
            key = lineL[0]
            assert key not in attrD, "Duplicate "+ key
            attrD[key] = {}
            for subkey, value in zip(attributeL, lineL[1:]):
                attrD[key][subkey] = value
    #-----------------------------------------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    header = 1
    aDict = {}
    for line in fh:
        lineL = line.strip().split('\t')
        if header:
            headerL = [i.upper() for i in lineL[1:]]
            for name in headerL:
                assert name not in aDict, "Duplicate "+name+'. STH WRONG!!!'
                aDict[name] = {}
            header -= 1
            continue
        #-----------------------------------------------
        key = lineL[0]
        for name, value in zip(headerL, lineL[1:]):
            assert key not in aDict[name], "Duplicate "+name+'. STH WRONG!!!'
            aDict[name][key] = value
    getSurvialData(aDict) #:--------------------------------------------------------------------------------
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    keyL = aDict.keys()
    #keyL.sort(key=lambda x: int(aDict[x]['Days.survival']))
    headerL = aDict[keyL[0]].keys()
    headerL.sort()
    #headerL.insert(0, 'ID')
    if attrD:
        print "%s\t%s\t%s\t%s\t%s" % ('ID', 'oldID', "SampleType", '\t'.join(attributeL), "\t".join(headerL))
    else:
        print "%s\t%s\t%s\t%s" % ('ID', 'oldID', "SampleType", "\t".join(headerL))
    for key in keyL:
        fullL = sampleIDd.get(key, [key])
        for each in fullL:
            if attrD:
                print "%s\t%s\t%s\t%s\t%s" % (each, key, typeD.get(each, 'unknown'), 
                        '\t'.join([attrD.get(each, {}).get(i, 'NA') for i in attributeL]), 
                        '\t'.join([aDict[key][i] for i in headerL]))
            else:
                print "%s\t%s\t%s\t%s" % (each, key, typeD.get(each, 'unknown'),  '\t'.join([aDict[key][i] for i in headerL]))

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


