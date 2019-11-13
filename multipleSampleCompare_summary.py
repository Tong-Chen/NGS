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
    This is designed to integrate the results of multipleSampleCompare.py.

    It will plot the number of DE genes for each compare, gene expression profiles of DE genes, pearson correlation of all samples, pca analysis of all samples, as well as annotated files for each group of comparing.

    All outputs are listed in a JSON file named **PREFIX.all.DE.summary.xls**.

'''

import sys
import os
#from json import dump as json_dump
import json
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import pandas as pd
#from multiprocessing.dummy import Pool as ThreadPool

debug = 0


def plot_de_count(all_DE, condL, type):
    count = all_DE+'.count.xls'
    count_fh = open(count, 'w')
    countD = {}
    for line in open(all_DE):
        type = line.split()[1]
        firstSep  = type.find('._')
        secondSep = type.find('_.')
        
        condA = type[:firstSep]
        cmp   = type[firstSep+2:secondSep]
        condB = type[secondSep+2:] 
        if cmp == 'higherThan':
            condA, condB = condB, condA
        else:
            assert cmp == 'lowerThan', cmp
        if condA not in countD:
            countD[condA] = {}
        if condB not in countD[condA]:
            countD[condA][condB] = 0
        countD[condA][condB] += 1
    #---------------------------------
    print >>count_fh, 'Samp\t'+'\t'.join(condL)
    for cond in condL:
        subD = countD.get(cond, {})
        tmpL = [str(subD.get(sec_cond, 'NA')) for sec_cond in condL]
        print >>count_fh, "{}\t{}".format(cond, '\t'.join(tmpL))
    count_fh.close()
    ## In this heatmap, each colored block represents 
    ## number of up-regulated DE genes or regions in samples in X-axis
    ## compared to samples in Y-axis.
    label = 'DE_{}_count'.format(type)
    cmd = ["s-plot heatmapS -f", count, "-A 45 -T 2 -l top -I",label,
          "-b TRUE -Y white"]
    if os.system(' '.join(cmd)):
        sys.exit(1)
#-----------------plot_de_count-------------
def plot_de_profile(all_DE, norm_mat_fl, condL, do_plot_de_profile):
    #print >>sys.stderr, norm_mat_fl
    #print >>sys.stderr, all_DE
    norm_mat = pd.read_table(norm_mat_fl, header=0, index_col=0)
    #print >>sys.stderr, norm_mat.iloc[0:3, 0:4]
    norm_mat.index.name = 'ID'
    norm_mat.index = norm_mat.index.map(str)
    #Extract DE_profile
    #all_DE_idL = list(set([line.split()[0] for line in open(all_DE)]))
    # 
    all_DE_idL = []
    all_DE_idL_top10 = []
    tmpD83 = {}
    for line in open(all_DE):
        gene, sample = line.split()
        all_DE_idL.append(gene)
        if sample not in tmpD83:
            tmpD83[sample] = 1
        if tmpD83[sample] < 11:
            all_DE_idL_top10.append(gene)
    all_DE_idL = list(set(all_DE_idL))
    all_DE_idL_top10 = list(set(all_DE_idL_top10))
    #print >>sys.stderr, all_DE_idL
    de_norm_mat = norm_mat[norm_mat.index.isin(all_DE_idL)]
    de_norm_mat_f = all_DE + '.norm' 
    de_norm_mat.index.name = 'ID'
    #print >>sys.stderr, de_norm_mat.iloc[0:3,0:4]
    de_norm_mat.to_csv(de_norm_mat_f, sep=b"\t")
    #Plot
    cmd = ['s-plot prettyHeatmap -f', de_norm_mat_f, '-c', str(len(condL))]
    if do_plot_de_profile:
        if os.system(' '.join(cmd)):
            sys.exit(1)

    de_norm_mat = norm_mat[norm_mat.index.isin(all_DE_idL_top10)]
    de_norm_mat_f = all_DE + '.norm.top10' 
    de_norm_mat.index.name = 'ID'
    #print >>sys.stderr, de_norm_mat.iloc[0:3,0:4]
    de_norm_mat.to_csv(de_norm_mat_f, sep=b"\t")
    #Plot
    cmd = ['s-plot prettyHeatmap -f', de_norm_mat_f, '-c', str(len(condL)), '-l \"-b TRUE\"']
    if do_plot_de_profile:
        if os.system(' '.join(cmd)):
            sys.exit(1)
    #os.system("/bin/rm -f "+de_norm_mat_f)
#-------------plot_de_profile-----------------

def plotPearson(norm_mat_fl, sampleFile):
    cmd = ['pearsonCorrelationMatrix.py -i', norm_mat_fl]
    os.system(' '.join(cmd))
    cmd = ["s-plot pheatmap -f", norm_mat_fl+".pearson.xls", 
           "-H TRUE -R TRUE -u 15 -v 15 -A 90", 
           "-t 'Pearson correlation of all samples'"]
    if sampleFile:
        cmd.extend(['-P', sampleFile, '-Q', sampleFile])
    print >>sys.stderr, ' '.join(cmd)
    os.system(' '.join(cmd))
#------------------------------------------------------
def pca(norm_mat_fl, sampleFile):
    cmd = ["s-plot pca -f", norm_mat_fl, "-L TRUE -w 22 -u 18"]
    if sampleFile:
        fh = open(sampleFile)
        header = fh.readline()
        Group = header.split()[1]
        fh.close()
        cmd.extend(["-g", sampleFile])
    else:
        Group = 'group'
    cmd.extend(["-c", Group, "-S", Group])
    print >>sys.stderr, ' '.join(cmd)
    os.system(' '.join(cmd))
#------------------------------------------


def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "\n\t%prog -c compare_pair -P prefix -a anno -n norm"
    parser = OP(usage=usages)
    parser.add_option("-c", "--compare-pair", dest="filein",
        help="compare_pair given to multipleSampleCompare.py.")
    parser.add_option("-s", "--sample-file", dest="sample_file",
        help="Sample file given to multipleSampleCompare.py. [Optional]")
    parser.add_option("-P", "--prefix", dest="prefix",
        help="[UPPERCASE P] The output prefix given to multipleSampleCompare.py.")
    parser.add_option("-n", "--norm-matrix", dest="norm_mat",
        help="Normalized data matrix file which will be used for \
sample classification,  sample correlation and DE profile analysis. \
Normally the file (if normalized) given to multipleSampleCompare.py \
or the count file given to multipleSampleCompare.py after normalization.")
    parser.add_option("-t", "--type", dest="type",
        help="Attributes of analyzed variables, currently <gene> or <region> are supported.")
    parser.add_option("-T", "--top", dest="top_regions",
        default=200000, type="int", 
        help="For results with hundreds of millions of regions, only display top <N> regions. Default=200000.")
    parser.add_option("-m", "--meltRegionBinMatrix", dest="meltRegionBinMatrix",
        help="If <type> is region, supplying parameters for <meltRegionBinMatrix.py> in format like <-b bin_size -c chromosome-file>. Unused now.")
    parser.add_option("-r", "--redo-DE", dest="DE",
        default=False, action="store_true", 
        help="Re-do DE analysis with given parameters.")
    parser.add_option("-p", "--no-plot-DE-profile", dest="no_plot_de_profile",
        default=False, action="store_true", 
        help="Specify to stop plotting DE profile. Normally should specify this parameter except very very many de regions.")
    parser.add_option("-C", "--no-cluster", dest="no_cluster_analysis",
        default=False, action="store_true", 
        help="Specify if sample cluster pca and pearson correlation should be done. Normally one should specify this parameter except very very many DE regions or part samples.")
    parser.add_option("-q", "--qvalue", dest="qvalue",
        default=0.1, type='float', help="FDR used for selecting DE genes. Default 0.1.If -r set to True, this will bed used to seclect DE regions or genes. If -r set to False, this value should be the same as qvalue used for DE analysis of input files.")
    parser.add_option("-Q", "--pvalue", dest="pvalue",
        default=0.05, type='float', help="P-value used for selecting DE genes. Default 0.05. If -r set to True, this will bed used to seclect DE regions or genes. If -r set to False, this value should be the same as pvalue used for DE analysis of input files.")
    parser.add_option("-f", "--log2fc", dest="log2fc",
        default=1, type='float', help="log2 fold change used for selecting DE genes. Default 1 meaning 2-fold change. If -r set to True,  this will bed used to seclect DE regions or genes. If -r set to False,  this value should be the same as log2fc used for DE analysis of input files.")
    parser.add_option("-a", "--anno", dest="anno",
        help="Annotation file. Optional. The first column of annotation should be unique and correspoding DE genes or regions name. One header line of annotation file is needed.")
    parser.add_option("-A", "--annoed", dest="annoed",
        default=False, action="store_true", help="Only used for output of <mergeMultipleSampleCompare.py> when all regions are annotated. Default false.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    assert options.type != None, "A string needed for -t"
    assert options.norm_mat != None, "A filename needed for -n"
    return (options, args)
#--------------------------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    compare_pair = options.filein
    sample_file  = options.sample_file
    prefix_g     = options.prefix
    DE_analysis  = options.DE
    verbose      = options.verbose
    fdr          = options.qvalue
    log2fc       = options.log2fc
    log2fc_neg   = (-1) * log2fc
    anno         = options.anno
    norm_mat_fl  = options.norm_mat
    type         = options.type
    annoed       = options.annoed
    pvalue       = options.pvalue
    top_regions  = options.top_regions
    do_cluster_analysis = not options.no_cluster_analysis
    do_plot_de_profile = not options.no_plot_de_profile
    up_count = 0
    dw_count = 0
    meltRegionBinMatrix = options.meltRegionBinMatrix
    global debug
    debug        = options.debug
    if anno:
        annoMat = pd.read_table(anno, header=0, index_col=0)
    #-----------------------------------
    condL = []
    if compare_pair == '-':
        fh = sys.stdin
    else:
        fh = open(compare_pair)
    #--------------------------------
    all_DE = prefix_g + '.all.DE'
    all_DE_summary = all_DE + '.summary.xls'
    all_DE_summaryL = []
    comp_pD = {}
    if DE_analysis:
        all_DE_fh = open(all_DE, 'w')
    for line in fh:
        condA, condB = line.strip().split()
        cond_key = '\t'.join([condA, condB])
        condL.append(condA)
        condL.append(condB)
        comp_pD[cond_key] = {}
        suffix = 'xls'
        if annoed:
            suffix = 'anno.xls'
        DE_results = '.'.join([prefix_g, condA+'._vs_.'+condB, suffix])
        if top_regions:
            top_regions_p1 = top_regions+1
            suffix_new = 'anno.top'+str(top_regions)+'.xls'
            DE_results_old = DE_results
            DE_results = '.'.join([prefix_g, condA+'._vs_.'+condB, suffix_new])
            cmd_top = ["head -n", str(top_regions_p1), 
                    DE_results_old, '>', DE_results]
            if os.system(' '.join(cmd_top)):
                print >>sys.stderr, " ".join(cmd_top)
                sys.exit(1)
        #----------------------------------------------
        comp_pD[cond_key]['all'] = DE_results
        if debug:
            print >>sys.stderr, DE_results
        if anno or DE_analysis:
            DE_results_matrix = pd.read_table(DE_results, header=0, index_col=0)
            DE_results_matrix.index.name = 'ID'
        if anno:
            DE_results_matrix_anno = DE_results_matrix.join(annoMat, how="left").sort_values('p_value')
            DE_results_anno = '.'.join([prefix_g, condA+'._vs_.'+condB, 'anno.xls'])
            DE_results_matrix_anno = DE_results_matrix_anno.fillna('-')
            DE_results_matrix_anno.index.name = 'ID'
            DE_results_matrix_anno.to_csv(DE_results_anno, sep=b'\t')
            comp_pD[cond_key]['all'] = DE_results_anno
        up_type = condA+'._lowerThan_.'+condB
        dw_type = condA+'._higherThan_.'+condB
        up = '.'.join([prefix_g, up_type, suffix])
        dw = '.'.join([prefix_g, dw_type, suffix])
        comp_pD[cond_key]['lowerThan'] = up
        comp_pD[cond_key]['higherThan'] = dw

        if DE_analysis:
            up_mat = DE_results_matrix[(DE_results_matrix['p_value']<=pvalue) & (DE_results_matrix['fdr']<=fdr) & (DE_results_matrix['log2FC']>=log2fc)]
            up_mat.to_csv(up, sep=b'\t')
            for i in up_mat.index.values:
                print >>all_DE_fh, "{}\t{}".format(i, up_type)
            dw_mat = DE_results_matrix[(DE_results_matrix['p_value']<=pvalue) & (DE_results_matrix['fdr']<=fdr) & (DE_results_matrix['log2FC']<=log2fc_neg)]
            dw_mat.to_csv(dw, sep=b'\t')
            for i in dw_mat.index.values:
                print >>all_DE_fh, "{}\t{}".format(i, dw_type)
        #----------------------------------------------
        #print >>sys.stderr, up
        if os.path.exists(up):
            up_matrix = pd.read_table(up, header=0, index_col=0)
            up_count = up_matrix.shape[0]
            #print >>sys.stderr, up_count
            up_matrix.index.name = 'ID'
        else:
            up_count = 0
            if anno:
                print >>sys.stderr, up+' not found'
                return 1
        if os.path.exists(dw):
            dw_matrix = pd.read_table(dw, header=0, index_col=0)
            dw_count = dw_matrix.shape[0]
            dw_matrix.index.name = 'ID'
        else:
            dw_count = 0
            if anno:
                print >>sys.stderr, dw+' not found'
                return 1
        comp_pD[cond_key]['lowerThan_Cnt'] = up_count
        comp_pD[cond_key]['higherThan_Cnt'] = dw_count
        if anno:
            up_matrix_anno = up_matrix.join(annoMat, how="left").sort_values('p_value')
            up_matrix_anno = up_matrix_anno.fillna('-')
            up_matrix_anno.index.name = 'ID'
            up_anno = '.'.join([prefix_g, up_type, 'anno.xls'])
            up_matrix_anno.to_csv(up_anno, sep=b'\t')
            dw_matrix_anno = dw_matrix.join(annoMat, how="left").sort_values('p_value')
            dw_matrix_anno = dw_matrix_anno.fillna('-')
            dw_matrix_anno.index.name = 'ID'
            dw_anno = '.'.join([prefix_g, dw_type, 'anno.xls'])
            dw_matrix_anno.to_csv(dw_anno, sep=b'\t')
            comp_pD[cond_key]['lowerThan'] = up_anno
            comp_pD[cond_key]['higherThan'] = dw_anno
        #-----------END one pair------------------------
    #-------END compare_pair-------------------------------
    if DE_analysis:
        all_DE_fh.close()
    #-------------------------------------
    condL = list(set(condL))
    condL.sort()
    plot_de_count(all_DE, condL, type)
    all_DE_summaryD = {}
    all_DE_summaryD['DE_count'] = {'file': all_DE+'.count.xls', 
            'plot': all_DE+'.count.xls.heatmapS.pdf'}
    #norm_mat = pd.read_table(norm_mat_fl, header=0, index_col=0)
    #norm_mat.index.name = 'ID'
    #norm_mat.index = norm_mat.index.map(str)
    #if do_plot_de_profile:
    plot_de_profile(all_DE, norm_mat_fl, condL, do_plot_de_profile)
    all_DE_summaryD['DE_profile'] = {'file': all_DE + '.norm.kmeans.xls', 
            'plot': all_DE + '.norm.kmeans.sort.heatmapS.pdf'}
    
    all_DE_summaryD['pearson'] = {'file': norm_mat_fl+".pearson.xls", 
        'plot':norm_mat_fl+".pearson.xls.pheatmap.pdf"}
    all_DE_summaryD['pca'] = {'plot':norm_mat_fl+".pca.scale.pdf"}
    if do_cluster_analysis:
        plotPearson(norm_mat_fl, sample_file)

        pca(norm_mat_fl, sample_file)
    #-----------------------------------------------

    all_DE_summaryD['DE_parameters'] = {'fdr': fdr, 'log2fc':log2fc, 'P_value': pvalue}

    all_DE_summaryL.append(all_DE_summaryD)
    all_DE_summaryL.append(comp_pD)
    
    with open(all_DE_summary, 'w') as all_DE_summary_fh:
        json.dump(all_DE_summaryL, all_DE_summary_fh, indent=4, sort_keys=True)
    #-------------END reading compare_pair----------
    #----close compare_pair handle for files-----
    if compare_pair != '-':
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
    ###---------procompare_pair the program---------
    #import procompare_pair
    #procompare_pair_output = sys.argv[0]+".prof.txt")
    #procompare_pair.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(procompare_pair_output)
    #p.sort_stats("time").print_stats()
    ###---------procompare_pair the program---------


