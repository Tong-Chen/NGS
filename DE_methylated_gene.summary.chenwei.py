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
    This is designed to summarize results output by DESeq2.sh and GO, KEGG enrichment file.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import re
from tools import *
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
reload(sys)
sys.setdefaultencoding('utf8')

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
    usages = "%prog -f file"
    parser = OP(usage=usages)
    parser.add_option("-f", "--files", dest="filein",
        metavar="FILEIN", help="The all.DE ")
    parser.add_option("-l", "--label", dest="label",
        help="One of F, M, FM")
    parser.add_option("-e", "--enrichmentType", dest="enrich_type",
        default="BP_GO, MF_GO, CC_GO, KEGG", 
        help="Default <BP_GO, MF_GO, CC_GO, KEGG>. And only these 4 types of enrichment are supported. These strings will be used to get filenames.")
    parser.add_option("-g", "--go-prefix", dest="go_prefix",
        help="String like <go/NK_trans.rc.xls.DESeq2.all.DE.entrez>.")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='3_expression_profile', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '3_expression_profile'.")
#    parser.add_option("-d", "--doc-only", dest="doc_only",
#        default=False, action="store_true", help="Specify to only generate doc.")
#    parser.add_option("-n", "--number", dest="number", type="int", 
#        default=40, help="Set the maximum allowed samples for barplot. Default 40.\
# If more than this number of samples are given, heatmap will be used instead. ")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def generateDoc(report_dir, report_sub_dir, file, go_prefix, enrich_typeL, label, curation_label):
    
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)

    print "\n### 自然生产和辅助生产小孩差异甲基化分析 {#dmc-%s}\n" % label
    
    knitr_read_txt(report_dir,  curation_label)
    

    prefix = file.replace('.all.DE', '')
    '''
    typeL = ['A._vs_.C_up', 'A._vs_.C_dw']
    '''
    typeL = list(set([line.split()[1] for line in open(file)]))
    '''
    typeD = {'A._vs_.C':['up', 'dw']}
    '''
    #typeD = {}
    #for cmp_type in typeL:
    #    cmp, type = cmp_type.split('_', 1)
    #    if cmp not in typeD:
    #        typeD[cmp] = [type]
    #    else:
    #        typeD[cmp].append(type)
    #typeL = typeD.keys()        
    #--------------------------------------------        
    #anno_all = prefix + '.normalized.anno.xls'
    #pca = prefix + ".normalized.rlog.pca.pdf"
    #pearson = prefix + ".normalized.rlog.pearson.pdf"
    #pca_png = prefix + ".normalized.rlog.pca.png"
    #pearson_png = prefix + ".normalized.rlog.pearson.png"
    Au_up_totalSum = prefix + ".Au_up.gene_diff_profile.totalSum.xls"
    Au_up_totalSum_png = prefix + ".Au_up.gene_diff_profile.totalSum.xls.pheatmap.png"
    Nu_up_totalSum = prefix + ".Nu_up.gene_diff_profile.totalSum.xls"
    Nu_up_totalSum_png = prefix + ".Nu_up.gene_diff_profile.totalSum.xls.pheatmap.png"

    copy(dest_dir, Au_up_totalSum, Nu_up_totalSum)
    copypng(dest_dir, Au_up_totalSum_png, Nu_up_totalSum_png)
    
    rAu_up_totalSum, rAu_up_totalSum_png, rNu_up_totalSum, rNu_up_totalSum_png = getRelativeDir([Au_up_totalSum, Au_up_totalSum_png, Nu_up_totalSum, Nu_up_totalSum_png], 
            report_sub_dir)

    print '''

(ref:dm-heatmap-{label}) 左图表示在辅助生殖后代 (**[Au]({rAu_up_totalSum})**)中启动子区甲基化水平增高的基因。右图表示在自然生产后代 (**[Nu]({rNu_up_totalSum})**) 中启动子区甲基化水平增高的基因。

```{{r dm-heatmap-{label}, out.width="50%", fig.cap="(ref:dm-heatmap-{label})"}}
knitr::include_graphics(c("{rAu_up_totalSum_png}", "{rNu_up_totalSum_png}"))
```
'''.format(label=label, rAu_up_totalSum=rAu_up_totalSum, rAu_up_totalSum_png=rAu_up_totalSum_png, rNu_up_totalSum=rNu_up_totalSum, rNu_up_totalSum_png=rNu_up_totalSum_png)
    
    print "\n### 差异修饰基因功能注释和富集分析"
    
    count = 0

    go_prefix_up = go_prefix+'.'+typeL[0]+'.'
    go_prefix_dw = go_prefix+'.'+typeL[1]+'.'
    enrich_cnt = 0
    for enrichType in enrich_typeL:
        enrich_cnt += 1
        go_prefix_up_enrich_xls = go_prefix_up+enrichType+'.xls'
        go_prefix_up_enrich_pdf = go_prefix_up+enrichType+'.scatterplot.dv.pdf'
        go_prefix_up_enrich_png = go_prefix_up+enrichType+'.scatterplot.dv.png'
        go_prefix_dw_enrich_xls = go_prefix_dw+enrichType+'.xls'
        go_prefix_dw_enrich_pdf = go_prefix_dw+enrichType+'.scatterplot.dv.pdf'
        go_prefix_dw_enrich_png = go_prefix_dw+enrichType+'.scatterplot.dv.png'
        copy(dest_dir, go_prefix_up_enrich_xls, go_prefix_dw_enrich_xls)
        copypdf(dest_dir, go_prefix_up_enrich_pdf, go_prefix_dw_enrich_pdf)

        rgo_prefix_up_enrich_xls, rgo_prefix_up_enrich_png, rgo_prefix_up_enrich_pdf, \
            rgo_prefix_dw_enrich_xls, rgo_prefix_dw_enrich_png, rgo_prefix_dw_enrich_pdf = \
            getRelativeDir([go_prefix_up_enrich_xls, go_prefix_up_enrich_png, go_prefix_up_enrich_pdf, \
            go_prefix_dw_enrich_xls, go_prefix_dw_enrich_png, go_prefix_dw_enrich_pdf], report_sub_dir)

        print '''
(ref:enrich-{label}-{count}-{enrich_cnt}) Top 30 gene ontology enrichment term or KEGG pathway for differentially methylated genes. Au represents assisted reproduction and Nu represents natural reproduction. 缺少的图表示对应样品中没有富集的条目。 Full lists of GO enrichment terms can be downloaded for [**{samp1}**]({rgo_prefix_up_enrich_xls}) ([ PDF pic](rgo_prefix_up_enrich_pdf)) and [**{samp2}**](rgo_prefix_dw_enrich_xls) ([ PDF pic](rgo_prefix_dw_enrich_pdf)).             

```{{r enrich-{label}-{count}-{enrich_cnt}, out.width="50%", fig.cap="(ref:enrich-{label}-{count}-{enrich_cnt})"}}
knitr::include_graphics(c("{rgo_prefix_up_enrich_png}", "{rgo_prefix_dw_enrich_png}"))
```

'''.format(label=label, count=count, enrich_cnt=enrich_cnt, samp1=typeL[0], samp2=typeL[1], rgo_prefix_up_enrich_xls=rgo_prefix_up_enrich_xls, rgo_prefix_up_enrich_pdf=rgo_prefix_up_enrich_pdf, rgo_prefix_dw_enrich_xls=rgo_prefix_dw_enrich_xls, rgo_prefix_dw_enrich_pdf=rgo_prefix_dw_enrich_pdf, rgo_prefix_up_enrich_png=rgo_prefix_up_enrich_png, rgo_prefix_dw_enrich_png=rgo_prefix_dw_enrich_png)
#-------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    label = options.label
    enrich_type = options.enrich_type
    enrich_typeL = re.split(r'[, ]*', enrich_type.strip())
    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    #doc_only = options.doc_only
    go_prefix = options.go_prefix
    curation_label = os.path.split(sys.argv[0])[-1].replace('.', '_')
    #if doc_only:
    #    generateDoc(report_dir, report_sub_dir, file, go_prefix, enrich_typeL, curation_label)
    #    return 0
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    aDict = {}
    generateDoc(report_dir, report_sub_dir, file, go_prefix, enrich_typeL, label, curation_label)

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


