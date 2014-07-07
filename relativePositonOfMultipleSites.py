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
    This is originally designed to check the relative position of a
    motif with miRNA target sites.

Motif_file: (GGACT part is used for filter only)

chr1    38076319        38076324        P_187_4@1  5.395544        +       GGACT
chr4    134972439       134972444       P_20432@1       4.571010    -       AAACT
chr4    134972345       134972350       P_20432@2       5.395544    -       GGACT
chr16   46448694        46448699        P_11118_1@1     4.607086    -       AGACT
chr16   46448671        46448676        P_11118_1@2     4.571010    -       AAACT
chr1    38076122        38076127        P_187_3@1       4.305981    +       TGACT
chr16   46448838        46448843        P_11118_2@1     4.607086    -       AGACT

Target_file: (mmu part is used for filter only)

chr17   56905154        56905175        P_13112_2:mmu-miR-29b-1-5p  158.00  +
chr16   3979445 3979466 P_10448_1:mmu-miR-29b-1-5p      164.00  -
chr16   3979423 3979445 P_10448_1:mmu-miR-29b-1-5p:2    156.00  -
chr7    25207631        25207652        P_25506_2:mmu-miR-29b-1-5p  158.00  +
chr4    106111466       106111490       P_19487_1:mmu-miR-29b-1-5p  158.00  +
chr13   52742453        52742473        P_7408_1:mmu-miR-29b-1-5p   161.00  +
chr12   16545241        16545264        P_6145_2:mmu-miR-29b-1-5p   166.00  -
chr5    93380012        93380033        P_21917_2:mmu-miR-29b-1-5p  174.00  +
chr10   76950769        76950792        P_2061_2:mmu-miR-29b-1-5p   156.00  -

'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i motif_file -t target_file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--motif-file", dest="m_f",
        metavar="MOTIF_FILE", help="A standard bed file. \
The sample file we used is the output of \
<transferfindMotifsOutputToBed.py>. We assume the regions \
in this file are smaller than those in the file given to -t.")
    parser.add_option("-s", "--m-separtor", dest="m_name_sep",
        default='@', help="The separtor for names in bed file. \
Default '@' for sample files.")
    parser.add_option("-f", "--m-filter", dest="m_filter_file",
        help="The files containing filters with each at one row. \
Only lines matched with filters will be used for following analysis.")
    parser.add_option("-t", "--target-file", dest="t_f",
        metavar="TARGET_FILE", help="A standard bed file. \
The sample file we used is the output of \
<transferMirandaOutputToBed.py>.")
    parser.add_option("-S", "--t-separtor", dest="t_name_sep",
        default=':', help="The separtor for names in bed file. \
Default ':' for sample files.")
    parser.add_option("-F", "--t-filter", dest="t_filter_file",
        help="The files containing filters with each at one row. \
Only lines matched with filters will be used for following analysis.")
    parser.add_option("-r", "--t-range", dest="t_range",
        help="Label specific regions for comparing, for example, \
seed regions of miRNAs (2-8) will be represented as 2-8 \
(1 started, both included, strand information will be considered \
innerly. Mltiple subregions can be given in format like \
2-8,10-14,15-17")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.m_f != None, "A filename needed for -i"
    assert options.t_f != None, "A filename needed for -t"
    return (options, args)
#--------------------------------------------------------------------

def computePos(m_start, m_end, t_start, t_end, t_strand, suffix):
    if m_end <= t_start:
        dist = t_start - m_end + 1
        if t_strand == '+':
            return "UP_"+str(dist) + suffix
        elif t_strand == '-':
            return "DW_"+str(dist) + suffix
    elif m_start >= t_end:
        dist = m_start - t_end + 1
        if t_strand == '+':
            return "DW_"+str(dist) + suffix
        elif t_strand == '-':
            return "UP_"+str(dist) + suffix
    else:
        if m_start <= t_start:
            overlap = m_end - t_start
            if t_strand == '+':
                return "OUP_"+str(overlap) + suffix
            elif t_strand == '-':
                return "ODW_"+str(overlap) + suffix
        elif m_end > t_end:
            overlap = t_end - m_start
            if t_strand == '+':
                return "ODW_"+str(overlap) + suffix
            elif t_strand == '-':
                return "OUP_"+str(overlap) + suffix
        else:
            return 'IN_1' + suffix
#---------------------------------------------------

def generateTargetPos(t_lineL, t_rangeL, revCompl=1):
    tmpD = {}
    tmpL = []
    t_start = int(t_lineL[1])
    t_end   = int(t_lineL[2])
    t_strand  = t_lineL[5]
    tmpD["TOTAL"] = [t_start, t_end, t_strand]
    tmpL.append("TOTAL")
    if t_rangeL:
        count = 0
        for subL in t_rangeL:
            count += 1
            sub_start = int(subL[0])
            sub_end   = int(subL[1])
            if (t_strand == '+' and revCompl) or \
                    (t_strand == '-' and (not revCompl)):
                new_end = t_end - sub_start + 1
                new_start = t_end - sub_end
            elif (t_strand == '-' and revCompl) or \
                    (t_strand == '+' and (not revCompl)):
                new_start = t_start + sub_start - 1
                new_end   = t_start + sub_end
            key = "SUB_" + str(count) 
            tmpD[key] = [new_start, new_end, t_strand]
            tmpL.append(key)
    #---------------------------------------------
    return tmpD, tmpL
#---------------------------------------------------

def relativePos(t_lineL, m_valueLL, t_rangeL):
    '''
    lineL = [chr, start, end, name.value, strand]
    m_valueLL = [ [chr1, start, end, name, value, strand, other], 
                  [bedL],
                  ...
                ]
    '''
    #compare whole region
    targetL, targetD = generateTargetPos(t_lineL, t_rangeL)
    #print m_valueLL
    for m_valueL in m_valueLL:
        #print m_valueL
        tmpL = t_lineL[:]
        tmpL.extend(m_valueL)
        m_start = int(m_valueL[1])
        m_end   = int(m_valueL[2])
        for target in targetL:
            t_start = targetL[target][0]
            t_end   = targetL[target][1]
            t_strand   = targetL[target][2]
            newtmp = computePos(m_start, m_end, t_start, t_end, t_strand,
                '@'+target)   
            tmpL.append(newtmp)
        #-----------------------------------
        print '\t'.join(tmpL)
#--------------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    m_file = options.m_f
    m_sep  = options.m_name_sep
    m_filter = options.m_filter_file
    t_file = options.t_f
    t_sep  = options.t_name_sep
    t_filter = options.t_filter_file
    t_range = options.t_range
    t_rangeL = []
    if t_range:
        for j in t_range.split(','):
            t_rangeL.append([int(i) for i in j.split('-')])
    #----------------------------------------
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    m_filter_l = []
    if m_filter:
        m_filter_l = [line.strip() for line in open(m_filter)]
    #--------------------------------
    m_dict = {}
    '''
    m_dict = {key:[[chr1, start, end, name, value, strand, other], 
                   [bedL],
                   ... 
                  ]
            }
    '''
    for line in open(m_file):
        lineL = line.split()
        if (not m_filter_l) or \
            (m_filter_l and lineL[-1] in m_filter_l):
            key = lineL[3].split(m_sep)[0]
            if key not in m_dict:
                m_dict[key] = [lineL]
            else:
                m_dict[key].append(lineL)
    #-------------END reading file----------
    t_filter_l = []
    if t_filter:
        t_filter_l = [line.strip() for line in open(t_filter)]
    #-----------------------------------------
    t_keyL = {}
    for line in open(t_file):
        lineL = line.split()
        key, filter = lineL[3].split(t_sep)[:2] 
        t_keyL[key] = 1
        if (not t_filter_l) or \
            (t_filter_l and filter in t_filter_l):
            if key in m_dict:
                relativePos(lineL, m_dict[key], t_rangeL)
            else:
                print '\t'.join(lineL), "No motif"
    #-------------END reading file----------
    for m_key, m_itemL in m_dict.items():
        if m_key not in t_keyL:
            for m_item in m_itemL:
                print '\t'.join(m_item), "Untargeted"
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



