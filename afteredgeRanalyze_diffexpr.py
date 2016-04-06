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
    This is designed to integrate the results of <DE_analysis.pl>,
    <edgeR.pl> and <analyze_diff_expr.pl>.
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
        metavar="FILEIN", help="compare_pair")
    parser.add_option("-D", "--dir", dest="dir",
        metavar="DIRECTORY", help="The directory containing all data")
    parser.add_option("-P", "--prefix", dest="prefix",
        help="The prefix of the project with <isoform> or <gene> \
specified, like $(prefix).isoform or $(prefix).gene.")
    parser.add_option("-p", "--pvalue", dest="pvalue",
        default=0.05, help="FDR used for selecting DE genes")
    parser.add_option("-f", "--log2fc", dest="log2fc",
        default=2, help="log2 fold change used for selecting DE genes")
    parser.add_option("-e", "--expr", dest="expr",
        help="<quantification/$(prefix).isoform.expr.TMM.fpkm.matrix>")
    parser.add_option("-a", "--anno", dest="anno",
        default="Trinotate/Trinotate_annotation_report.xls", 
        help="<Trinotate/Trinotate_annotation_report.xls>")
    parser.add_option("-I", "--annoIndex", dest="index",
        default=2, help="For <isoforms>, index should be 2 \
representing indexing annotation by the second column ids. \
For <genes>, index should be 1 meaning indexing annotation \
by the first column.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def readEnriched(*enriched):
    header = 1
    for label, file in enriched:
        file_out = file + '.xls'
        label += "-UP.enriched"
        out_fh = open(file_out, 'w')
        header = 1
        for line in open(file):
            lineL = line.strip().split('\t')
            newLineL = [lineL[8], lineL[1], lineL[3], lineL[7]]
            if header:
                newLineL.append("Sample")
                header -= 1
            else:
                newLineL.append(label)
            print >>out_fh, '\t'.join(newLineL)
        out_fh.close()
    #--------------------------------------
#--------------------------------------

def readDepleted(*depleted):
    header = 1
    for label, file in depleted:
        file_out = file + '.xls'
        label += "-UP.depleted"
        out_fh = open(file_out, 'w')
        header = 1
        for line in open(file):
            lineL = line.strip().split('\t')
            newLineL = [lineL[9], lineL[2], lineL[3], lineL[8]]
            if header:
                newLineL.append("Sample")
                header -= 1
            else:
                newLineL.append(label)
            print >>out_fh, '\t'.join(newLineL)
        out_fh.close()
    #--------------------------------------
#--------------------------------------

def getTopTerm(type, prefix, top, *file):
    file_out_n = prefix+".DE_genes.GOseq."+type+'.xls'
    file_out = open(file_out_n, 'w')
    print >>file_out, "Term\tneg_log10pvalue\tCount\tFDR\tSample"
    for single in file:
        count = 0
        for line in open(single):
            if line.startswith(type):
                print >>file_out, line,
                count += 1
            if count >= top:
                break
    #---------------------------------------
    file_out.close()
    cmd = ['s-plot scatterplotDoubleVariable -f', file_out_n, 
        '-o Sample -v Term -c neg_log10pvalue -s Count -w 40 -a 70 -E pdf -R 90 -H 0 -V 1']
    os.system(' '.join(cmd))
    convert = ['convert -density 150 -quality 90',
            file_out_n+'.scatterplot.dv.pdf',
            file_out_n+'.scatterplot.dv.png']
    os.system(' '.join(convert))
#------------------------------------

def readAnno(anno, index):
    header = 1
    annoD = {}
    for line in open(anno):
        line = line.strip()
        if header:
            head = line
            header -= 1
        else:
            lineL = line.split('\t', 3)
            annoD[lineL[index]] = line
    #--------------------------
    return annoD, head
#-----------------------------
def readExpr(expr):
    header = 1
    exprD = {}
    for line in open(expr):
        line = line.strip()
        if header:
            head = 'expr_TMM_FPKM'+line
            header -= 1
        else:
            lineL = line.split('\t', 2)
            exprD[lineL[0]] = line
    #----------------
    return exprD, head
#---------------------------
def annoDE_results(DE_results, annoD, annoH, exprD, exprH):
    file_out = DE_results + '.anno.xls'
    if debug:
        print >>sys.stderr, file_out
    file_fh = open(file_out, 'w')
    header = 1
    for line in open(DE_results):
        line = line.strip()
        if header:
            print >>file_fh, '\t'.join(["ID", line, exprH, annoH])
            header -= 1
        else:
            key = line.split('\t', 1)[0]
            #print >>sys.stderr, key
            #print >>sys.stderr, exprD.keys()
            #print >>sys.stderr, annoD.keys()

            print >>file_fh, '\t'.join([line, exprD[key],
                annoD.get(key, "No annotation")])
    file_fh.close()
#_----------------------------------------

def annoSubset(annoD, annoH, *subset):
    for file, label in subset:
        file_out =  file.replace('subset', 'xls')
        if debug:
            print >>sys.stderr, file_out
        file_fh = open(file_out, 'w')
        header = 1
        for line in open(file):
            line = line.strip()
            if header:
                print >>file_fh, '\t'.join([line, annoH])
                header -= 1
            else:
                key = line.split('\t', 1)[0]
                print >>file_fh, '\t'.join([line, 
                    annoD.get(key, "No annotation")])
                print '%s\t%s' % (key, label.replace('-', '_'))
        file_fh.close()
#_----------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    compare_pair = options.filein
    dir = options.dir
    verbose = options.verbose
    pvalue = options.pvalue
    log2fc = options.log2fc
    prefix_g = options.prefix
    global debug
    debug = options.debug
    anno = options.anno
    index = int(options.index) - 1
    annoD, annoH = readAnno(anno, index)
    expr = options.expr
    exprD, exprH = readExpr(expr)
    #-----------------------------------
    if compare_pair == '-':
        fh = sys.stdin
    else:
        fh = open(compare_pair)
    #--------------------------------
    for line in fh:
        condA, condB = line.strip().split()
        DE_results = '.'.join([dir+'/'+prefix_g+'.expr.counts.matrix', 
            condA+'_vs_'+condB, 'edgeR.DE_results'])
        if debug:
            print >>sys.stderr, DE_results
        annoDE_results(DE_results, annoD, annoH, exprD, exprH)

        MA = DE_results + 'MA.png'
        volcano = DE_results + 'Volcano.png'
        prefix = '.'.join([DE_results, 'P'+pvalue+'_C'+log2fc])
        subsetA = '.'.join([prefix, condA+'-UP', 'subset'])
        subsetB = '.'.join([prefix, condB+'-UP', 'subset'])
        
        labelA = condA+'_vs_'+condB+'.'+condA+'_UP'
        labelB = condA+'_vs_'+condB+'.'+condB+'_UP'

        annoSubset(annoD, annoH, [subsetA, labelA], [subsetB, labelB])

        subsetA_enriched = subsetA + '.GOseq.enriched'
        subsetA_depleted = subsetA + '.GOseq.depleted'

        subsetB_enriched = subsetB + '.GOseq.enriched'
        subsetB_depleted = subsetB + '.GOseq.depleted'

        readEnriched([labelA, subsetA_enriched], [labelB,
            subsetB_enriched])
        #readDepleted([labelA, subsetA_depleted], [labelB,
        #    subsetB_depleted])
        
        getTopTerm('BP', prefix, 20, subsetA_enriched+'.xls',
            subsetB_enriched+'.xls')

        getTopTerm('MF', prefix, 20, subsetA_enriched+'.xls',
            subsetB_enriched+'.xls')

        getTopTerm('CC', prefix, 20, subsetA_enriched+'.xls',
            subsetB_enriched+'.xls')

#        getTopTerm('BP', prefix, 50, subsetA_enriched+'.xls',
#            subsetA_depleted+'.xls', subsetB_enriched+'.xls',
#            subsetB_depleted+'.xls')
#
#        getTopTerm('MF', prefix, 50, subsetA_enriched+'.xls',
#            subsetA_depleted+'.xls', subsetB_enriched+'.xls',
#            subsetB_depleted+'.xls')
#
#        getTopTerm('CC', prefix, 50, subsetA_enriched+'.xls',
#            subsetA_depleted+'.xls', subsetB_enriched+'.xls',
#            subsetB_depleted+'.xls')
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


