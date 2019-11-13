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
from tools import *
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
    parser.add_option("-D", "--dir", dest="dir",
        metavar="DIR", help="The directory containing all data")
    parser.add_option("-T", "--target-dir", dest="report_dir",
        metavar="TARGET-DIR", help="The directory in <report> for saving data.")
    parser.add_option("-P", "--prefix", dest="prefix_pro",
        metavar="PREFIX", help="The prefix of the project with gene or \
isoforms specified. LIke Andro.gene or Andro.isoform")
    parser.add_option("-n", "--num_de_genes", dest="num_de_genes",
        help="<numDE_feature_counts>. The program will \
generate <numDE_feature_counts.Ppvalue_Clof2fc.matrix.heatmapS.png>.")
    parser.add_option("-H", "--heatmap", dest="heatmap",
        metavar="FILEIN",
        help="<$(prefix).matrix.trend.4.kmeans>")
    parser.add_option("-L", "--title", dest="title",
        help="A string to name this document.")
    parser.add_option("-t", "--type", dest="type",
        metavar="TYPE", help="Accept `Unigenes` or `isoforms`")
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="compare_pair")
    parser.add_option("-p", "--pvalue", dest="pvalue",
        default=0.05, help="FDR used for selecting DE genes")
    parser.add_option("-f", "--log2fc", dest="log2fc",
        default=2, help="log2 fold change used for selecting DE genes")
    parser.add_option("-k", "--keggEnrichment", dest="kegg_enrichment_prefix",
        help="Files given to goEnrichment_long.py with tags added like <Saussurea_involucrata.DEG.KEGG>.")
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
        '-o Sample -v Term -c neg_log10pvalue -s Count -w 20 -a 20 -E pdf -R 90 -H 0 -V 1']
    os.system(' '.join(cmd))
    convert = ['convert -density 150 -quality 90',
            file_out_n+'.scatterplot.dv.pdf',
            file_out_n+'.scatterplot.dv.png']
    os.system(' '.join(convert))
#------------------------------------



def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    compare_pair = options.filein
    dir = options.dir
    report_dir = options.report_dir + '/'
    prefix_pro = options.prefix_pro
    type = options.type
    title = options.title
    verbose = options.verbose
    pvalue = options.pvalue
    log2fc = options.log2fc
    kegg_enrichment_prefix = options.kegg_enrichment_prefix

    num_de_genes_dat = dir + '/' +options.num_de_genes + '.P' + pvalue + "_C" \
        + log2fc + '.matrix'
    num_de_genes = num_de_genes_dat + '.heatmapS.pdf'
    heatmap_dat = dir + '/' + options.heatmap + '/' + \
        options.heatmap + ".result"
    heatmap = heatmap_dat  + ".final.sort.xls.heatmapS.pdf"
    global debug
    debug = options.debug
    #-----------------------------------
    os.system('mkdir -p report/'+report_dir)
    copy('report/'+report_dir, num_de_genes_dat, heatmap_dat)
    copypdf('report/'+report_dir, num_de_genes)
    copypdf('report/'+report_dir, heatmap)
    print "# {} \n\n".format(title)
    
    curation_label = "DE_{}_summary_profile".format(type)
    knitr_read_txt('report', curation_label)

    print "## Number of differentially expressed (DE) %s\n" % type

    print """
(ref:de-num-gene) Heatmap showing number of differentially expressed %s at log2 fold change(log2FC) larger than **%s** and false discovery rate (FDR) less than **%s**. 差异表达基因数目统计。 灰颜色代表对应样品没有做差异基因比较分析。[PDF](%s) [SOURCE](%s)

```{r de-num-gene, fig.cap="(ref:de-num-gene)"}
knitr::include_graphics("%s")
```
""" % (type, log2fc, pvalue, report_dir+os.path.split(num_de_genes)[1], 
         report_dir+os.path.split(num_de_genes_dat)[1], 
         report_dir+os.path.split(num_de_genes)[1].replace('pdf', 'png'))

    print "## Expression profile of differentially expressed (DE) %s\n" % type

    print """
(ref:expr-profile) Summary patterns of differentially expressed %s at log2 fold change(log2FC) larger than **%s** and false discovery rate (FDR) less than **%s**. Each line represents on gene or transcript. 差异基因表达谱热图展示。[PDF](%s) [SOURCE](%s)

```{r expr-profile, fig.cap="(ref:expr-profile)"}
knitr::include_graphics("%s")
```
"""  % (type, log2fc, pvalue, report_dir+os.path.split(heatmap)[1],
        report_dir+os.path.split(heatmap_dat)[1], 
        report_dir+os.path.split(heatmap)[1].replace('pdf', 'png'))

    #-----------------------------------
    if compare_pair == '-':
        fh = sys.stdin
    else:
        fh = open(compare_pair)
    #--------------------------------
    print "\n## 样品两两比较差异基因分析\n"

    curation_label = "DE_gene_summary_each_two_sample"
    knitr_read_txt('report', curation_label)

    count = 0
    for line in fh:
        count += 1
        condA, condB = line.strip().split()
        DE_results = '.'.join([dir+'/'+prefix_pro+'.counts.matrix', 
            condA+'_vs_'+condB, 'edgeR.DE_results'])
        #MA = DE_results + '.MA.png'
        #volcano = DE_results + '.Volcano.png'
        MA_volcano = DE_results + '.MA_n_Volcano.pdf'
        prefix = '.'.join([DE_results, 'P'+pvalue+'_C'+log2fc])
        subsetA = '.'.join([prefix, condA+'-UP', 'subset'])
        subsetB = '.'.join([prefix, condB+'-UP', 'subset'])

        annoxls = DE_results + '.anno.xls'
        subsetAxls = '.'.join([prefix, condA+'-UP', 'xls'])
        subsetBxls = '.'.join([prefix, condB+'-UP', 'xls'])

        subsetA_enriched = subsetA + '.GOseq.enriched'
        subsetA_depleted = subsetA + '.GOseq.depleted'

        subsetB_enriched = subsetB + '.GOseq.enriched'
        subsetB_depleted = subsetB + '.GOseq.depleted'
        
        print "### Summary of DE %s between _%s_ and _%s_\n" % (type, condA, condB)
        #copypng('report/'+report_dir, MA, volcano)
        copypdf('report/'+report_dir, MA_volcano)
        MA_png = MA_volcano.replace('.pdf', '-0.png')
        volcano_png = MA_volcano.replace('.pdf', '-1.png')
        copypng('report/'+report_dir, MA_png, volcano_png)

        print """
(ref:MA-volcano-DE-%d) **Left**. MA plot of differentially expressed (DE) %s \
at log2 fold change (log2FC) larger than %s and false discovery rate (FDR) less than %s. MA-plots are used to study dependences between the log ratio of two variables and the mean values of the same two variables. The log ratios of the two measurements are called M values (from 'minus' in the log scale) and are represented in the vertical axis. The mean values of the two measurements are called A values (from 'average' in the log scale) and are represented in the horizontal axis. MA plot用来查看表达数值与差异倍数分布，其Y轴中值应位于0点处。红色的点代表差异显著的基因或转录本。**Right** Volcano plot of differentially expressed (DE) %s at log2 fold change (log2FC) larger than **%s** and false discovery rate (FDR) less than **%s**. A volcano plot is a scatter plot that is often used when giving an overview of interesting genes or regions. The log fold change is plotted on the x-axis and the negative log10 p-value is plotted on the y-axis. Positive logFC represents %s up-regulated in %s. Negative logFC represents %s up-regulated in %s. 火山图用来查看差异基因的分布。X轴的负半轴表示在样品%s中上调，正半轴表示在样品%s中上调。红色的点代表差异显著的基因或转录本。上下图展示的是同一个意思，下面的图的点更小些。[pdf](%s)

```{r MA-volcano-DE-%d, out.width="48%%", fig.cap="(ref:MA-volcano-DE-%d)"}
knitr::include_graphics(c("%s", "%s"))
```

""" % (count, type, log2fc, pvalue, type, log2fc, pvalue, type, 
        condB, type, condA, condA, condB, 
        report_dir+os.path.split(MA_volcano)[1], 
        count, count, 
        report_dir+os.path.split(MA_png)[1], 
        report_dir+os.path.split(volcano_png)[1])
        #report_dir+os.path.split(MA)[1], 
        #report_dir+os.path.split(volcano)[1])

        copy('report/'+report_dir, annoxls, subsetAxls, subsetBxls)
        print "### Annotation of DE \
%s between _%s_ and _%s_\n" % (type, condA, condB)
        #print "[ct-columns]\n"
        #print "[ct-column=0.9]\n"
        print "\n**Annotation files 差异基因注释结果:**\n"
        print "* [%s_vs_%s.anno.xls](%s/%s)\n" % \
            (condA, condB, report_dir, os.path.split(annoxls)[1])
        print "* [%s-UP.anno.xls](%s/%s)\n" % \
            (condA, report_dir, os.path.split(subsetAxls)[1])
        print "* [%s-UP.anno.xls](%s/%s)\n" % \
            (condB, report_dir, os.path.split(subsetBxls)[1])

        print "\n\n**GO Enrichment files 差异基因功能富集结果:**\n"
        
        for  cond in [condA, condB]:
            #for type in ['enriched', 'depleted']:
            for go_type in ['enriched']:
                tmpFile = '.'.join([prefix, cond+'-UP',
                    'subset.GOseq', go_type, 'xls'])
                copy('report/'+report_dir, tmpFile)
                print "* [%s-UP.subset.GOseq.%s.xls](%s/%s)\n" % \
                    (cond, go_type, report_dir, os.path.split(tmpFile)[1])
        #------------------------------
        #print "[/ct-columns]\n"
        print 

        if kegg_enrichment_prefix:
            kegg_enrichment_prefix2 = kegg_enrichment_prefix + '.' + condA+'_vs_'+condB
            print "\n\n**KEGG Enrichment files 差异基因功能富集结果:**\n"
            for  cond in [condA, condB]:
                #for type in ['enriched', 'depleted']:
                for go_type in ['enriched']:
                    tmpFile = '.'.join([kegg_enrichment_prefix2, 
                        cond+'_UP.xls'])
                    copy('report/'+report_dir, tmpFile)
                    print "* [%s-UP.KEGG.%s.xls](%s/%s)\n" % \
                        (cond, go_type, report_dir, os.path.split(tmpFile)[1])
            #------------------------------
            #print "[/ct-columns]\n"
            print 

        BP_dat = prefix + ".DE_genes.GOseq.BP.xls"
        BP = BP_dat + ".scatterplot.dv.pdf"
        CC_dat = prefix + ".DE_genes.GOseq.CC.xls"
        CC = CC_dat + ".scatterplot.dv.pdf"
        MF_dat = prefix + ".DE_genes.GOseq.MF.xls"
        MF = MF_dat + ".scatterplot.dv.pdf"

        #CC = prefix + ".DE_genes.GOseq.CC.xls.scatterplot.dv.png"
        #MF = prefix + ".DE_genes.GOseq.MF.xls.scatterplot.dv.png"
        copypdf('report/'+report_dir, BP, CC, MF)
        copy('report/'+report_dir, BP_dat, CC_dat, MF_dat)

        print "#### GO enrichment (BP) of DE %s between _%s_ and _%s_\n" \
                % (type, condA, condB)
        print """
(ref:GO-enrichment-BP%d) Enriched biological process. `neg_log10pvalue` represents the enrichment significance for each GO term. It is computed using formula `(-1)*(log10(multiple-testing corrected pvalue))`. `Count` represents the number of DE genes or transcripts enriched in each GO term. [PDF](%s) [SOURCE](%s)

```{r GO-enrichment-BP%d, out.width="99%s", fig.cap="(ref:GO-enrichment-BP%d)"}
knitr::include_graphics("%s")
```
""" % (count, report_dir+os.path.split(BP)[1],
               report_dir+os.path.split(BP_dat)[1], 
               count, '%', count, 
               report_dir+os.path.split(BP)[1].replace('pdf', 'png'))

        print "#### GO enrichment (MF) of DE %s between _%s_ and _%s_\n" \
                % (type, condA, condB)
        print """
(ref:GO-enrichment-MF%d) Enriched molecular functions. `neg_log10pvalue` represents the enrichment significance for each GO term. It is computed using formula `(-1)*(log10(multiple-testing corrected pvalue))`. `Count` represents the number of DE genes or transcripts enriched in each GO term. [PDF](%s) [SOURCE](%s)

```{r GO-enrichment-MF%d, out.width="99%s", fig.cap="(ref:GO-enrichment-MF%d)"}
knitr::include_graphics("%s")
```
""" % (count, report_dir+os.path.split(MF)[1],
               report_dir+os.path.split(MF_dat)[1], 
               count, '%', count, 
               report_dir+os.path.split(MF)[1].replace('pdf', 'png'))

        print "#### GO enrichment (CC) of DE %s between _%s_ and _%s_\n" \
                % (type, condA, condB)
        print """
(ref:GO-enrichment-CC%d) Enriched cellular component. `neg_log10pvalue` represents the enrichment significance for each GO term. It is computed using formula `(-1)*(log10(multiple-testing corrected pvalue))`. `Count` represents the number of DE genes or transcripts enriched in each GO term. [PDF](%s) [SOURCE](%s)

```{r GO-enrichment-CC%d, out.width="99%s", fig.cap="(ref:GO-enrichment-CC%d)"}
knitr::include_graphics("%s")
```
""" % (count, report_dir+os.path.split(CC)[1],
               report_dir+os.path.split(CC_dat)[1], 
               count, '%', count, 
               report_dir+os.path.split(CC)[1].replace('pdf', 'png'))

        if kegg_enrichment_prefix:
            kegg_dat1 = kegg_enrichment_prefix2 + '.' + condA + '_UP.xls'
            kegg_dat2 = kegg_enrichment_prefix2 + '.' + condB + '_UP.xls'
            kegg_dat1_pdf = kegg_dat1 + '.top.scatterplot.dv.pdf' 
            kegg_dat2_pdf = kegg_dat2 + '.top.scatterplot.dv.pdf' 
            copypdf('report/'+report_dir, kegg_dat1_pdf, kegg_dat2_pdf)
            copy('report/'+report_dir, kegg_dat1, kegg_dat2)

            print "#### KEGG enrichment of DE %s between _%s_ and _%s_\n" \
                    % (type, condA, condB)
            print """
(ref:GO-enrichment-KEGG%d) Enriched KEGG pathway. `neg_log10pvalue` represents the enrichment significance. It is computed using formula `(-1)*(log10(multiple-testing corrected pvalue))`. `Count` represents the number of DE genes or transcripts enriched in each term. [PDF1](%s) [PDF2](%s) [SOURCE1](%s) [SOURCE2](%s)

```{r GO-enrichment-KEGG%d, out.width="48%s", fig.cap="(ref:GO-enrichment-KEGG%d)"}
knitr::include_graphics(c("%s", "%s"))
```
""" % (count, report_dir+os.path.split(kegg_dat1_pdf)[1],
               report_dir+os.path.split(kegg_dat2_pdf)[1],
               report_dir+os.path.split(kegg_dat1)[1], 
               report_dir+os.path.split(kegg_dat2)[1], 
               count, '%', count, 
               report_dir+os.path.split(kegg_dat1_pdf)[1].replace('pdf', 'png'),
               report_dir+os.path.split(kegg_dat2_pdf)[1].replace('pdf', 'png'))
    #----------Clustering information------------------------
    #heatmap = dir + '/' + options.heatmap + '/' + \
    #    options.heatmap + ".result.final.sort.heatmapS.png"
#        readEnriched([condA, subsetA_enriched], [condB,
#            subsetB_enriched])
#        readDepleted([condA, subsetA_depleted], [condB,
#            subsetB_depleted])
#        
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


