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
    usages = "%prog -i \
quantification/trinity_quantification.DE/Sau.matrix.trend.16.kmeans \
-c 16 -p 0.05 -f 2"
    parser = OP(usage=usages)
    parser.add_option("-D", "--dir", dest="dir",
        metavar="DIR", help="The directory containing all data. \
Normally \
<quantification/trinity_quantification.DE/$(prefix).matrix.trend.16.kmeans>.")
    parser.add_option("-T", "--tar-dir", dest="tar_dir",
        metavar="TAR-DIR", help="The directory in <doc> for saving data.")
    parser.add_option("-c", "--cluster-num", dest="cluster_num",
        metavar="CLUSTER-NUM", help="The number of clusters.")
    parser.add_option("-p", "--pvalue", dest="pvalue",
        default=0.05, help="FDR used for selecting DE genes")
    parser.add_option("-f", "--log2fc", dest="log2fc",
        default=2, help="log2 fold change used for selecting DE genes")
    parser.add_option("-L", "--title", dest="title",
        metavar="TYPE", help="A string to represent this doc.")
    parser.add_option("-t", "--type", dest="type",
        metavar="TYPE", help="Accept `Unigenes` or `isoforms`")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.dir != None, "A filename needed for -D"
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



def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    dir = options.dir.rstrip('/')
    type = options.type
    title = options.title
    tar_dir = options.tar_dir.rstrip('/') + '/'
    prefix = dir.split('/')[-1]
    dir += '/'
    cluster_num = int(options.cluster_num) + 1
    verbose = options.verbose
    pvalue  = options.pvalue
    log2fc  = options.log2fc
    global debug
    debug = options.debug
    #-----------------------------------
    target_dir = "report/"+tar_dir
    os.system('mkdir -p '+target_dir)
    #copypdf(target_dir, num_de_genes)
    #copypng(target_dir, heatmap)
    smooth_sta  = prefix + ".cluster.mean.xls"
    smooth_line = smooth_sta + ".lines.smooth.pdf"
    copypdf(target_dir, dir+smooth_line)
    copy(target_dir, dir+smooth_sta)
    
    print '# {}\n'.format(title)
    
    curation_label = "Co_expr_patterns_{}".format(type) 
    knitr_read_txt("report", curation_label)

    print '''
## Expression pattern of DE {type}

(ref:expr-pat-{type}) Expression patterns of differentially expressed (DE) {type} at log2 fold change(log2FC) larger than **{log2fc}** and false discovery rate (FDR) less than **{pvalue}**. [PDF]({pdf}) [SOURCE]({xls})

```{{r expr-pat-{type}, fig.cap="(ref:expr-pat-{type})"}}
knitr::include_graphics("{png}")
```
'''.format(type=type, log2fc=log2fc, pvalue=pvalue, pdf=tar_dir+smooth_line,
           xls=tar_dir+smooth_sta, png=tar_dir+smooth_line.replace('pdf','png'))

    for i in range(1, cluster_num):
        i = str(i)
        line_dat = '.'.join([prefix, 'cluster', i, 'lines.xls'])
        line_plot = '.'.join([line_dat, 'lines.smooth.pdf'])
        copypdf(target_dir, dir+line_plot)
        copy(target_dir, dir+line_dat)
        print '## Expression pattern and functional enrichment \
analysis for %s in cluster _%s_\n' % (type, i)
        print '''
(ref:cluster-{cluster}-{type}) Expression pattern for cluster **{cluster}**. [PDF]({pdf}) [SOURCE]({xls})        

```{{r cluster-{cluster}-{type}, fig.cap="(ref:cluster-{cluster}-{type})" }}
knitr::include_graphics("{png}")
```
'''.format(cluster=i, pdf=tar_dir+line_plot, xls=tar_dir+line_dat, 
               png=tar_dir+line_plot.replace('pdf','png'), type=type)

        annoxls = '.'.join([prefix, i, 'anno.xls'])
        copy(target_dir, dir+annoxls)
        print "### annotation files\n"
        print "* [%s](%s/%s)\n" % (annoxls, tar_dir, annoxls)

        for go_type in ['GOseq.enriched.xls']:
            tmpFile = '.'.join([prefix, i, go_type])
            copy(target_dir, dir+tmpFile)
            print "* [%s](%s/%s)\n" % (tmpFile, tar_dir, tmpFile)

        #-----------------------------------
        BP_dat = '.'.join([prefix, i, 'GOseq.BP.xls'])
        BP = '.'.join([BP_dat, 'scatterplot.dv.pdf'])
        MF_dat = '.'.join([prefix, i, 'GOseq.MF.xls'])
        MF = '.'.join([MF_dat, 'scatterplot.dv.pdf'])
        CC_dat = '.'.join([prefix, i, 'GOseq.CC.xls'])
        CC = '.'.join([CC_dat, 'scatterplot.dv.pdf'])
        #MF = '.'.join([prefix, i, 'GOseq.MF.xls.scatterplot.dv.pdf'])
        #CC = '.'.join([prefix, i, 'GOseq.CC.xls.scatterplot.dv.pdf'])
        copypdf(target_dir, dir+BP, dir+MF, dir+CC)
        copy(target_dir, dir+BP_dat, dir+MF_dat, dir+CC_dat)
        print "### GO enrichemnt (BP) of DE %s in cluster _%s_\n" % (type, i)
        print '''
(ref:cluster-{cluster}-{type}-BP) Enriched biological process. `neg_log10pvalue` represents the enrichment significance for each GO term. It is computed using formula `(-1)*(log10(multiple-testing corrected pvalue))`. `Count` represents the number of DE genes or transcripts enriched in each GO term. [PDF]({pdf}) [SOURCE]({xls})

```{{r cluster-{cluster}-{type}-BP, fig.cap="(ref:cluster-{cluster}-{type}-BP)" }}
knitr::include_graphics("{png}")
```
'''.format(cluster=i, type=type, pdf=tar_dir+BP, xls=tar_dir+BP_dat, 
           png=tar_dir+BP.replace('pdf', 'png'))

        print "### GO enrichemnt (MF) of DE %s in cluster _%s_\n" % (type, i)

        print '''
(ref:cluster-{cluster}-{type}-MF) Enriched molecular function. `neg_log10pvalue` represents the enrichment significance for each GO term. It is computed using formula `(-1)*(log10(multiple-testing corrected pvalue))`. `Count` represents the number of DE genes or transcripts enriched in each GO term. [PDF]({pdf}) [SOURCE]({xls})

```{{r cluster-{cluster}-{type}-MF, fig.cap="(ref:cluster-{cluster}-{type}-MF)" }}
knitr::include_graphics("{png}")
```
'''.format(cluster=i, type=type, pdf=tar_dir+MF, xls=tar_dir+MF_dat, 
           png=tar_dir+MF.replace('pdf', 'png'))

        print "### GO enrichemnt (CC) of DE %s in cluster _%s_\n" % (type, i)

        print '''
(ref:cluster-{cluster}-{type}-CC) Enriched cellular component. `neg_log10pvalue` represents the enrichment significance for each GO term. It is computed using formula `(-1)*(log10(multiple-testing corrected pvalue))`. `Count` represents the number of DE genes or transcripts enriched in each GO term. [PDF]({pdf}) [SOURCE]({xls})

```{{r cluster-{cluster}-{type}-CC, fig.cap="(ref:cluster-{cluster}-{type}-CC)" }}
knitr::include_graphics("{png}")
```
'''.format(cluster=i, type=type, pdf=tar_dir+CC, xls=tar_dir+CC_dat, 
           png=tar_dir+CC.replace('pdf', 'png'))

        print 
    #----------Clustering information------------------------
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


