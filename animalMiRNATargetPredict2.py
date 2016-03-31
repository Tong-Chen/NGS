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
'''


import sys
import os
import string
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP



def readtarget(target):
    targetDict = {}
    for line in open(target):
        if line[0] == '>':
            locus = line[1:-1]
            assert locus not in targetDict
            targetDict[locus] = ''
        else:
            targetDict[locus] += line.strip().upper().replace('T','U')
    #-------------------------------------------
    return targetDict
#-------------------------------------------------

def alignOneMirnaMultipleTarget(targetDict, seed, edit_d, lv):
    #if edit_d == 0:
    #    target_pre = []
    #    mismatched = []
    #    for key, seq in targetDict.items():
    #        if seq.count(seed) >= hitNum:
    #            target_pre.append(key)    
    #    #----------------------------------------
    #    return target_pre, mismatched
    ##--------------------------------------
    #else:
    target_pre = {}
    mismatched = {}
    lenseed = len(seed)
    tmpseed = 'N' * lenseed
    for key, seq in targetDict.items():
        if 0:
            print '-----------------------'
        pos40 = seq.find(seed)
        mismatchedCnt = 0
        seedCount = seq.count(seed)
        if 0:
            print 'seq', seq
            print 'seed', seed
            print 'seedCount', seedCount
        if seedCount:
            seq2 = seq.replace(seed,tmpseed)
            if 0:
                print 'seq2', seq2
            seq = seq2
        #--------------------------------------
        if edit_d > 0:
            lenseq = len(seq)
            jump = 0
            for i in range(lenseq):
                start = i
                if jump:
                    jump -= 1
                    continue
                end = i + lenseed
                if (end>lenseq+1):
                    break
                elif end > lenseq:
                    end = lenseq
                tmpRegion = seq[start:end]
                if lv.distance(tmpRegion,seed) <= edit_d:
                    if 0:
                        print 'tmpRegion', tmpRegion
                        print 'seed', seed
                    mismatchedCnt += 1
                    #mismatched.append()
                    jump = lenseed
            #--------------------------
        #---------------------------
        target_pre[key] = seedCount
        mismatched[key] = mismatchedCnt
    #----------------------------------------
    return target_pre, mismatched
#----------------------------------------

def readMir(mirna):
    #------------------------------
    #>cel-miR-1-5p MIMAT0020301 Caenorhabditis elegans miR-1-5p
    #CAUACUUCCUUACAUGCCCAUA
    #------------------------------
    mirnaDict = {}
    for line in open(mirna):
        if line[0] == '>':
            locus = line[1:].strip().split()[0]
            mirnaDict[locus] = ''
        else:
            mirnaDict[locus] += line.strip().upper()
    #-------------------------------
    return mirnaDict
#--------------------------------------------
    



def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    desc = "Predict miRNA target (Animal). Output result to file(s)"
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-m", "--miRNA", dest="mir",
        metavar="mir.fa", help="A fasta format file. Usually \
downloaded from miRbase. The name of the miRNA is the part \
of string before the first blank with '>' removed.")
    parser.add_option("-t", "--targetSeq", dest="target",
        metavar="target.fa", help="The sequence file in fasta format \
used to be the potential target sequences.")
    parser.add_option("-e", "--mis-match", dest="edit_d",
        metavar=0, default=0, help="The max allowed mismatch. Default 0 means \
perfect match. If a number larger than 0 are given, perfect match will \
be first searched and then masked by N to search for mismatched-target \
site, ")
    parser.add_option("-n", "--num-of-seed", dest="seed_n",
        metavar=1,default=1, help="Only taken as target if no less than \
<seed_n> matched regions are found. Default 1.")
    parser.add_option("-s", "--start_pos", dest="start_pos",
        default=2, help="The start of seed regions in miRNA (1-based). \
Default 2.")
    parser.add_option("-l", "--end_pos", dest="end_pos",
        default=8, help="The end of seed regions in miRNA (1-based). \
Default 8.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.mir != None, "A filename needed for -m"
    assert options.target != None, "A filename needed for -t"
    return (options, args)
#--------------------------------------------------------------------


def main():
    '''
    1.character case, all transfer to uppercase
    2.all U and T will transfer to A in mirna seed reverse
    complementary.
    3.all utr sequence T->U
    '''
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    mirna = options.mir
    target = options.target
    edit_d = int(options.edit_d)
    if edit_d > 0:
        import Levenshtein as lv
    else:
        lv = ''
    start_pos = int(options.start_pos)-1
    end_pos = int(options.end_pos)
    hitNum = int(options.seed_n)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    transTable = string.maketrans('ACGUT', 'UGCAA')
    targetDict = readtarget(target)
    mirnaDict = readMir(mirna)
    #---mirna index----------------------
    #newmirnaDict = {} #for target index
    file = sys.argv[4]+ '.' + str(hitNum) + '.miRNA.index.target' 
    fh = open(file, 'w')
    for key, value in mirnaDict.items():
        value = value[start_pos:end_pos]
        value = value.translate(transTable)
        value = value[::-1]
        #newmirnaDict[key] = value
        target_pre, mismatched = alignOneMirnaMultipleTarget\
            (targetDict, value,edit_d, lv)
        keyL169 = target_pre.keys()
        #print target_pre
        for key169 in keyL169:
            per_match_n = target_pre[key169]
            mis_match_n = mismatched[key169]
            if per_match_n + mis_match_n >= hitNum:
                print >>fh, '%s\t%s\t%s\t%d\t%d' % \
                    (key, key169, value, per_match_n, mis_match_n)
    fh.close()
    #-------targetIndex--------------------------
    #file = sys.argv[4] + '.target.index.mirna' 
    #fh = open(file, 'w')
    #for key, value in targetDict.items():
    #    mirna_pre = []
    #    for mirna_key, mirna_value in newmirnaDict.items():
    #        count88 = value.count(mirna_value)
    #        if count88:
    #            mirna_pre.append(mirna_key+'\t'+str(count88))
    #        #-------------------------------------
    #    #-------------------------------------
    #    for item in mirna_pre:
    #        print >>fh, "%s\t%s" % (key, item)
    #fh.close()






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



