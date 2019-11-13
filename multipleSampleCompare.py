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
    This is designed to compare multiple samples with or without replications.


data_matrix
ID  A   B   C   D   E   F   G   H
a   1   2   1   2   1   2   1   2
b   1   2   1   2   1   2   1   2
c   1   2   1   2   1   2   1   2
d   1   2   1   2   1   2   1   2

sampleFile (Normally, two columns would be OK. If nested group needed, please add extra columns.)

Samp    GrpA    GrpB
A   Grp1    Cond1
B   Grp2    Cond2
C   Grp1    Cond1
D   Grp1    Cond3
E   Grp2    Cond2
F   Grp3    Cond3
G   Grp4    Cond1
H   Grp4    Cond2

compare_pair
Grp1    Grp2
Grp3    Grp4
Cond1   Cond2
Cond1   Cond3

norm_factor
A   1000
B   1000
C   1000
D   1000
E   1000
F   1000
G   1000
H   1000

fisher:
    YES

t.test:
  calculates the T-test for the means of TWO INDEPENDENT samples of scores.

wilcox.test (wilcoxon):
    The Wilcoxon signed-rank test tests the null hypothesis that two related paired samples come from the same distribution. In particular,  it tests whether the distribution of the differences x - y is symmetric about zero. It is a non-parametric version of the paired T-test. Because the normal approximation is used for the calculations,  the samples used should be large. A typical rule is to require that n > 20.

kw.test (kruskal):
    The Kruskal-Wallis H-test tests the null hypothesis that the population median of all of the groups are equal. It is a non-parametric version of ANOVA. The test works on 2 or more independent samples,  which may have different sizes. Note that rejecting the null hypothesis does not indicate which of the groups differs. Post-hoc comparisons between groups are required to determine which groups are different.

ks.test:
    Computes the Kolmogorov-Smirnov statistic on 2 samples.
    This is a two-sided test for the null hypothesis that 2 independent samples are drawn from the same continuous distribution.

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from math import log
from statsmodels.stats.multitest import multipletests
from scipy import stats

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
    usages = "\n%prog -i file\nzcat gzfile | %prog -i -"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A data matrix. Gzipped file should be zcat to STDIN.")
    parser.add_option("-s", "--sample-file", dest="sampleFile",
        metavar="Samp-cond", help="Format as described above to specify group information for each sample. If no sampleFile specified, each sample will be treated as one group.")
    parser.add_option("-c", "--compare-file", dest="compareFile",
        metavar="compare-pair", help="Format as described above to specify which two samples needed to be compared. If no compareFile specified, each two groups will be compared.")
    parser.add_option("-m", "--statistical-method", dest="method",
        default='t.test', 
        help="1. If each of the compared group contains only one sample, <fisher.test> will be used for count data. No statistical method will be used for float data.\n\n2. If one of compared groups contains more than one samples, <t.test> will be used in default. <wilcox.test, kw.test, ks.test> are also acceptable.")
    parser.add_option("-C", "--count-data", dest="count_data",
        default = False, action="store_true", 
        help="Uppercase C. Only effective when both groups containing only one sample. Default False. If fisher.test wanted, this should be set to True.")
    parser.add_option("-n", "--norm-factor", dest="norm_factor",
        help="A two columns file as described above to specify normalize factor for each sample. No normalization will be performed if not supplied.")
    parser.add_option("-S", "--scale-factor", dest="scale_factor",
        default=1, type='int', help="An integer to scale data to avoid too small values. Default 1 meanin no scale.")
    parser.add_option("-L", "--log2-transformed", dest="log2_already",
        default=False, action="store_true", 
        help="Specify if input file data has been already log2 transformed.. Default False.")
    parser.add_option("-l", "--log2-transform", dest="log2",
        default=False, action="store_true", 
        help="log2 transform data before comparing. Default False.")
    parser.add_option("-d", "--data-type", dest="data_type",
        help="A short string to describe the quantification data for output label, such as <count>, <CPM (when -n and -S is given), <FPKM>, <RPKM> or others. If <-L> is specified, the program will automatically add <Log2> like <Log2 CPM>.")
    parser.add_option("-o", "--output-prefix", dest="op_prefix",
        help="Prefix of output file.")
    parser.add_option("-F", "--log2FC", dest="log2FC",
        default=1, type="float", help="Log2 fold change (default 1, meaning 2-fold difference)")
    parser.add_option("-p", "--fdr", dest="fdr",
        default=0.1, type="float", help="Adjusted p-value (default 0.1)")
    parser.add_option("-O", "--no-override", dest="ovride",
        default=True, action="store_false", help="Force override already existed comparasion (Default True, when set to False the program will check if same compare has been done. Any changes in original matrix, norm file, compared pair will induce override.). If set to True, all comparision will be redone (which is the DEFAULT). Currently this should be always set to TRUE.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the programe")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    assert options.data_type != None, "A description string needed for -d"
    return (options, args)
#--------------------------------------------------------------------
def readSampleFile(sampleFile, sampleL):
    '''
    Samp    Conds   Grps
    A   Grp1    Cond1
    B   Grp2    Cond2
    C   Grp1    Cond3
    D   Grp1    Cond3
    E   Grp2    Cond2
    F   Grp3    Cond2
    G   Grp4    Cond1
    H   Grp4    Cond3
    '''
    sampleD = {}
    if sampleFile:
        header = 1
        for line in open(sampleFile):
            if header:
                header -= 1
                continue
            lineL = line.split()
            samp  = lineL[0]
            grpL  = lineL[1:]
            for grp in grpL:
                if grp not in sampleD:
                    sampleD[grp] = [samp]
                else:
                    sampleD[grp].append(samp)
                #---------------
            #-----------------
        #------------------
    else:
        sampleD = dict([i, [i]] for i in sampleL)
    return sampleD
#---------------------------------------------
def readCom_pair(compareFile, grpL):
    '''
    Grp1    Grp2
    Grp3    Grp4
    '''
    if compareFile:
        return [tuple(line.split()) for line in open(compareFile)]
    else:
        grpL.sort()
        len_grpL = len(grpL)
        compareL = []
        for i in range(len_grpL-1):
            for j in range(i+1, len_grpL):
                compareL.append((grpL[i], grpL[j]))
    return compareL
#---------------------------------------------
def read_norm_factor(norm_factor, sampleL):
    if norm_factor:
        return dict([[line.split()[0], float(line.split()[1])] 
                     for line in open(norm_factor)])
    else:
        return dict([[i, 1] for i in sampleL])

#------------------------------
def fisher(valueD, total_first, total_second, scale):
    '''
    valueD = {id1:[[1, 2], [3, 4]], id2:[[1, 3], [5, 6]]}
    '''
    warn = 1
    from fisher import pvalue
    correctL = []
    pList = []
    no_correctL = []
    for id, valueL in valueD.items():
        if debug:
            print id, valueL
        try:
            q = int(valueL[0][0])
            m = int(valueL[1][0])
        except ValueError:
            if warn:
                print >>sys.stderr, "Force transfer float number to integer!!!"
                warn = 0
            q = int(float(valueL[0][0]))
            m = int(float(valueL[1][0]))
        if q * m == 0:
            q = q+1
            m = m+1
        p = pvalue(q, m, total_first-q, total_second-m)
        p = p.two_tail
        q = q*scale/total_first
        m = m*scale/total_second
        tmpL = [id, q, m]
        diff = log(m / q, 2)
        tmpL.append(diff)
        tmpL.append(p)
        if abs(diff)>=0.2 and p < 0.35: 
            tmpL.append(1)
            correctL.append(tmpL)
            pList.append(p)
        else:
            tmpL.append(1)
            no_correctL.append(tmpL)
    if debug:
        print >>sys.stderr, pList
    if pList:
        p_adjL = multipletests(pList, method="fdr_bh")[1]
        for tmpL, p_adj in zip(correctL, p_adjL):
            tmpL[-1] = p_adj
    resultL = correctL[:]
    resultL.extend(no_correctL)
    resultL.sort(key=lambda x: x[-1])
    return resultL
#---------------------------------------------------
def direct_compare(valueD, total_first, total_second, scale):
    resultL = []
    for id, valueL in valueD.items():
        q = float(valueL[0][0])
        m = float(valueL[1][0])
        if q * m == 0:
            q = q+1
            m = m+1
        tmpL = [id, q, m]
        q = q*scale/total_first
        m = m*scale/total_second
        diff = log(m / q, 2)
        tmpL.append(diff)
        resultL.append(tmpL)
    resultL.sort(key=lambda x: x[-1])
    return resultL
#-----------------------------------------------------
def stat_pvalue(v1, v2, method):
    if method == 't.test':
        return stats.ttest_ind(v1, v2, equal_var=False)[1]
    elif method == 'wilcox.test':
        return stats.wilcoxon(v1, v2)[1]
    elif method == 'kw.test':
        return stats.kruskal(v1, v2)[1]
    elif method == 'ks.test':
        return stats.ks_2samp(v1, v2)[1]
#-----------------------------------------------------

def rep_compare(valueD, total_first, total_second, method, log2, scale, log2_already):
    resultL = []
    correctL = []
    pList = []
    no_correctL = []
    for id, valueL in valueD.items():
        tmpL = [id]
        #meanL = [id]
        if log2:
            v1 = [log((float(i)+1)*scale/total_first, 2) for i in valueL[0]]
        else:
            v1 = [float(i)*scale/total_first for i in valueL[0]]
        tmpL.extend(v1)
        len_v1 = len(v1)
        meanV1 = sum(v1) / len_v1
        if log2:
            v2 = [log((float(i)+1)*scale/total_second, 2) for i in valueL[1]] 
        else:
            v2 = [float(i)*scale/total_second for i in valueL[1]] 
        len_v2 = len(v2)
        meanV2 = sum(v2) / len_v2
        tmpL.extend(v2)
        tmpL.append(meanV1)
        tmpL.append(meanV2)
        #meanL.append()
        if meanV1 * meanV2 == 0:
            meanV1 += 1
            meanV2 += 1
        if log2_already:
            diff = meanV2-meanV1
        else:
            diff = log(meanV2/meanV1, 2)
        #if log2:
        #    v1 = [log(i+1, 2) for i in v1]
        #    v2 = [log(i+1, 2) for i in v2]
        tmpL.append(diff)
        if abs(diff) >= 0.2:
            p = stat_pvalue(v1, v2, method)
        else:
            p = 0.5
        tmpL.append(p)
        if abs(diff)>=0.2 and p < 0.2: 
            tmpL.append(1)
            correctL.append(tmpL)
            pList.append(p)
        else:
            tmpL.append(1)
            #no_correctL.append(tmpL)
    if pList:
        p_adjL = multipletests(pList, method="fdr_bh")[1]
        for tmpL, p_adj in zip(correctL, p_adjL):
            tmpL[-1] = p_adj
    resultL = correctL[:]
    #resultL.extend(no_correctL)
    resultL.sort(key=lambda x: x[-1])
    return resultL
    
#-----------------------------------------------------
def output_result(headerL, resultL, len_grp1, len_grp2, name, op_prefix, op_all_fh, log2FC, fdr):
    samp1, samp2 = name
    all_tag = '._vs_.'.join(name)
    all = op_prefix+'.'+all_tag+'.xls'
    all_fh = open(all, 'w')
    all_mean = op_prefix+'.'+all_tag+'.mean.xls'
    all_mean_fh = open(all_mean, 'w')
    sampDE = op_prefix+'.'+all_tag+'.all.DE'
    sampDE_fh = open(sampDE, 'w')
    up = '._higherThan_.'.join(name)
    up_fh = open(op_prefix+'.'+up+'.xls', 'w') 
    dw = '._lowerThan_.'.join(name)
    dw_fh = open(op_prefix+'.'+dw+'.xls', 'w') 
    print >>all_fh, '\t'.join(headerL)
    if len_grp1*len_grp2<=1:
        print >>all_mean_fh, '\t'.join(lineL)
    else:
        headerL2 = [headerL[0]]
        headerL2.extend(headerL[len_grp1+len_grp2+1:])
        print >>all_mean_fh, '\t'.join(headerL2)
    #---------------------------------------------
    print >>up_fh, '\t'.join(headerL)
    print >>dw_fh, '\t'.join(headerL)
    diff = headerL.index('log2FC')
    len_headerL = len(headerL)
    if debug:
        print >>sys.stderr, headerL
        print >>sys.stderr, diff
    for lineL in resultL:
        if debug:
            print >>sys.stderr, lineL
        tmp_log2FC = lineL[diff]
        if len_headerL == diff+1:
            padj = 0
        else:
            padj = lineL[-1]
        lineL = [str(i) for i in lineL]
        print >>all_fh, '\t'.join(lineL)
        if len_grp1*len_grp2<=1:
            print >>all_mean_fh, '\t'.join(lineL)
        else:
            lineL2 = [lineL[0]]
            lineL2.extend(lineL[len_grp1+len_grp2+1:])
            print >>all_mean_fh, '\t'.join(lineL2)
        #---------------------------------------------
        if padj < fdr:
            if tmp_log2FC >= log2FC:
                print >>dw_fh, '\t'.join(lineL)
                print >>op_all_fh, '\t'.join([lineL[0], dw])
                print >>sampDE_fh, '\t'.join([lineL[0], dw])
            elif tmp_log2FC <= (-1)*log2FC:
                print >>up_fh, '\t'.join(lineL)
                print >>op_all_fh, '\t'.join([lineL[0], up])
                print >>sampDE_fh, '\t'.join([lineL[0], up])
    #------------------------------------------------
    all_fh.close()
    all_mean_fh.close()
    up_fh.close()
    dw_fh.close()
    sampDE_fh.close()
#-------------------------------------------------------

def ovrideF(ovride, name, op_prefix, headerL):
    if ovride:
        return 1
    else:
        all = '._vs_.'.join(name)
        all = op_prefix+'.'+all+'.xls'
        if os.path.exists(all):
            all_fh = open(all)
            header = all_fh.readline().strip()
            all_fh.close()
            if '\t'.join(headerL) == header:
                return 0
        return 1
#----------------------------------

def DE_analysis(compareD, sampD, normD, count_data, scale, log2, 
        op_prefix, method, op_all_fh, log2FC, fdr, ovride, log2_already, data_type):
    '''
    compareD = {(grp1, grp2):{id1:[[1, 2], [3, 4]], id2:[[1, 3], [5, 6]]}}
    '''
    for compare_pair, valueD in compareD.items():
        name = compare_pair
        grp1, grp2 = compare_pair
        headerL = ['ID']
        sampL1 = sampD[grp1]
        if debug:
            print >>sys.stderr, sampL1
        len_grp1 = len(sampL1)
        sampL2 = sampD[grp2]
        if debug:
            print >>sys.stderr, sampL2
        len_grp2 = len(sampL2)
        headerL.extend([tmp_samp+'('+data_type+')' for tmp_samp in sampL1])
        headerL.extend([tmp_samp+'('+data_type+')' for tmp_samp in sampL2])
        if len_grp1*len_grp2 == 1:
            if count_data:
                headerL.extend(['log2FC', 'p_value', 'fdr'])
                ovride_tmp = ovrideF(ovride, name, op_prefix, headerL)
                if not ovride_tmp:
                    continue
                resultL = fisher(valueD, normD[sampL1[0]], normD[sampL2[0]], scale)
            else:
                headerL.append('log2FC')
                ovride_tmp = ovrideF(ovride, name, op_prefix, headerL, scale)
                if not ovride_tmp:
                    continue
                resultL = direct_compare(valueD, normD[sampL1[0]], normD[sampL2[0]])
        else:
            headerL.extend([grp1+'(mean)', grp2+'(mean)', 'log2FC', 'p_value', 'fdr'])
            ovride_tmp = ovrideF(ovride, name, op_prefix, headerL)
            if not ovride_tmp:
                continue
            resultL = rep_compare(valueD,normD[sampL1[0]],normD[sampL2[0]],method,log2,scale, log2_already)
        #-------------------------------------------------------------------
        valueD = None
        output_result(headerL, resultL, len_grp1, len_grp2, name, op_prefix, op_all_fh, log2FC, fdr)
        #compareD.pop(compare_pair)
        resultL = None

#------------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file         = options.filein
    sampleFile   = options.sampleFile
    compareFile  = options.compareFile
    count_data   = options.count_data
    norm_factor  = options.norm_factor
    scale_factor = options.scale_factor
    log2         = options.log2
    log2_already = options.log2_already
    op_prefix    = options.op_prefix
    log2FC       = options.log2FC
    fdr          = options.fdr
    op_all       = op_prefix+'.all.DE'
    op_all_fh    = open(op_all, 'w')
    verbose      = options.verbose
    ovride       = options.ovride
    method       = options.method
    data_type    = options.data_type
    
    if log2_already:
        log2 = False

    if log2:
        log2_already = True

    if log2 or log2_already:
        data_type = "Log2 "+data_type

    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    header = 1
    samp2posD = {}
    '''
    sampleD = {'grp1':[A, B], 'grp2':[C, D]}
    compareL = [(grp1, grp2), (grp1, grp3)]
    normD = {'A':1, 'B':1, 'C':1}
    compareD = {(grp1, grp2):{id1:[[1, 2], [3, 4]], id2:[[1, 3], [5, 6]]}}
    '''
    compareD = {}
    for line in fh:
        lineL = line.strip().split('\t')
        if header:
            header -= 1
            len_lineL = len(lineL)
            for i in range(1, len_lineL):
                samp2posD[lineL[i]] = i
                if debug:
                    print >>sys.stderr, lineL[i], i
            sampleD  = readSampleFile(sampleFile, samp2posD.keys())
            compareL = readCom_pair(compareFile, sampleD.keys())
            for compare_pair in compareL:
                compareD[compare_pair] = {}
            normD    = read_norm_factor(norm_factor, samp2posD.keys())
            #print samp2posD["THP1_DMSO_15"]
            if debug:
                print >>sys.stderr, "samp2posD", samp2posD
                print >>sys.stderr, "sampleD", sampleD
                print >>sys.stderr, "compareL", compareL
                print >>sys.stderr, "normD", normD
                print >>sys.stderr, "compareD", compareD

            continue
        #------------------------------------------------------
        id = lineL[0]
        for compare_pair in compareL:
            grp1, grp2 = compare_pair
            value1 = [lineL[samp2posD[samp]] for samp in sampleD[grp1]]
            value2 = [lineL[samp2posD[samp]] for samp in sampleD[grp2]]
            compareD[compare_pair][id] = [value1, value2]
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    DE_analysis(compareD, sampleD, normD, count_data, scale_factor, log2, 
        op_prefix, method, op_all_fh, log2FC, fdr, ovride, log2_already, data_type)
    op_all_fh.close()
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


