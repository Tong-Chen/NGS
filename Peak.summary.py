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
    This is designed to summarize peak calling results.

    Result files are combined with <samp_name> <type> <gene_attr>

        Peak_file: <samp_name>_<type>/<samp_name>_<type>.peak.bed
        BW file: <samp_name>_<type>/<samp_name>_<type>.log2ratio.bw
        
        EnrichXLS: <summary_dir>/<type>.<gene_attr>.all.entrez.<samp_name>_<type>.<enrich_type>.scatterplot.dv.xls
        EnrichPDF: <summary_dir>/<type>.<gene_attr>.all.entrez.<samp_name>_<type>.<enrich_type>.scatterplot.dv.pdf
        EnrichPNG: <summary_dir>/<type>.<gene_attr>.all.entrez.<samp_name>_<type>.<enrich_type>.scatterplot.dv.png
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
    parser.add_option("-s", "--sample-name", dest="samp_name",
        metavar="FILEIN", help="The name of samples. Normally $(sampleL).")
    parser.add_option("-t", "--type", dest="type",
        help="Different types of modification. Normally $(typeL).")
    parser.add_option("-g", "--gene-attr", dest="gene_attr",
        default="peak_gene_promoter",
        help="Different ways of finding target genes. Normally support <peak_gene_promoter>, <peak_gene_closest>.")
    parser.add_option("-e", "--enrichmentType", dest="enrich_type",
        default="BP_GO, MF_GO, CC_GO, KEGG", 
        help="Default <BP_GO, MF_GO, CC_GO, KEGG>. And only these 4 types of enrichment are supported. These strings will be used to get filenames.")
    parser.add_option("-d", "--go-dir", dest="go_dir",
        help="Directory containing GO enrichment results. Normally $(summary_dir).")
    parser.add_option("-T", "--tag", dest="tag",
        default='tag', help="A string containing only aphabets to use as program specific label. Normally this is unused unless this program is used multiple times in same pipeline. Default <tag>.")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='5_peak_profile', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '5_peak_profile'.")
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
    assert options.samp_name != None, "A filename needed for -s"
    return (options, args)
#--------------------------------------------------------------------

def generateDoc(report_dir, report_sub_dir, tag, go_dir,
        sample_nameL, typeL, gene_attrL, enrich_typeL):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)

    print "\n# Peak鉴定和修饰图谱 {#peak-identification-modi-profile}\n"
    
    knitr_read_txt(report_dir, "Peak_identification_profile")
    
    count = 0
    for type in typeL:
        print "\n## " + type + ' profile\n\n'
        for samp in sample_nameL:
            count += 1
            print "\n### {} {}\n\n".format(samp, type)
            full_name = samp+'_'+type
            peak = full_name+'/'+full_name+'.peak.bed'
            bw   = full_name+'/'+full_name+'.log2ratio.bw'
            copy(dest_dir, peak, bw)
            peak,bw = getRelativeDir([peak,bw], report_sub_dir)
            peak_n, bw_n = getFileName([peak, bw])
            print "* Peak file in bed format: [{}]({})\n".format(peak_n, peak)
            #print "* Peak bound or modified gene list in text format: [{}]({})\n".format(gene_n, gene)
            print "* Track file in bigWig format: [{}]({})\n".format(bw_n, bw)

            if 'peak_gene_promoter' in gene_attrL:
                middle = "peak_gene_promoter.all.entrez"
                tag2 = tag + 'peak-gene-promoter'
                print "\n### 启动子区存在修饰的基因\n"
                gene = full_name+'/'+full_name+'.peak_gene_promoter.list'
                copy(dest_dir, gene)
                print "\n* 启动子区存在修饰的基因: [XLS]({})".format(getRelativeDir(gene, report_sub_dir))
                print "\n#### 样品{}{}修饰或结合基因(启动子区)的功能富集分析\n".format(samp, type)

                print '''
{} 修饰或结合的基因的筛选标准是：修饰或结合区域落在基因的启动子区域（转录起始位点上游1kb、下游 500 nt)。

'''.format(type)
        
                go_prefix = '.'.join([go_dir+'/'+type, middle, full_name+'.'])
                enrich_cnt = 0
                for enrichType in enrich_typeL:
                    enrich_cnt += 1
                    go_prefix_enrich_xls = go_prefix+enrichType+'.xls'
                    go_prefix_enrich_pdf = go_prefix+enrichType+'.scatterplot.dv.pdf'
                    copy(dest_dir, go_prefix_enrich_xls)
                    copypdf(dest_dir, go_prefix_enrich_pdf)
                    go_prefix_enrich_xls,go_prefix_enrich_pdf = \
                        getRelativeDir([go_prefix_enrich_xls,go_prefix_enrich_pdf],
                                report_sub_dir)
                    go_prefix_enrich_xls_name,go_prefix_enrich_pdf_name = \
                        getFileName([go_prefix_enrich_xls,go_prefix_enrich_pdf])

                    print '''
(ref:enrich-{tag}-{count}-{enrich_cnt}) Top 30 gene ontology or KEGG enrichment terms for genes regulated by {type} in {samp}. 缺少的图表示对应样品中没有富集的条目。 Full lists of GO or KEGG enrichment terms can be downloaded ([{go_prefix_enrich_xls_name}]({go_prefix_enrich_xls})) ([{go_prefix_enrich_pdf_name}]({go_prefix_enrich_pdf})).             

```{{r enrich-{tag}-{count}-{enrich_cnt}, out.width="90%", fig.cap="(ref:enrich-{tag}-{count}-{enrich_cnt})"}}
knitr::include_graphics(c("{go_prefix_enrich_png}"))
```

'''.format(tag=tag2, count=count, enrich_cnt=enrich_cnt, type=type, samp=samp, 
        go_prefix_enrich_xls_name=go_prefix_enrich_xls_name, 
        go_prefix_enrich_pdf_name=go_prefix_enrich_pdf_name, 
        go_prefix_enrich_xls=go_prefix_enrich_xls, 
        go_prefix_enrich_pdf=go_prefix_enrich_pdf, 
        go_prefix_enrich_png=go_prefix_enrich_pdf.replace('pdf', 'png'))

            if 'peak_gene_closest' in gene_attrL:
                middle = "peak_gene_closest.all.entrez"
                tag2 = tag + 'peak-gene-closest'
                print "\n### 修饰区邻近的基因\n"
                gene = full_name+'/'+full_name+'.peak_gene_closest.list'
                copy(dest_dir, gene)
                print "\n* 修饰区临近基因列表: [XLS]({})\n".format(getRelativeDir(gene, report_sub_dir))
                print "\n#### 样品{}{}修饰或结合基因(启动子区)的功能富集分析\n".format(samp, type)

                print '''
{} 修饰区临近基因的筛选标准是：修饰或结合区域最近的基因 (上游优先, 距离记载在最后一列)。

'''.format(type)
        
                go_prefix = '.'.join([go_dir+'/'+type, middle, full_name+'.'])
                enrich_cnt = 0
                for enrichType in enrich_typeL:
                    enrich_cnt += 1
                    go_prefix_enrich_xls = go_prefix+enrichType+'.xls'
                    go_prefix_enrich_pdf = go_prefix+enrichType+'.scatterplot.dv.pdf'
                    copy(dest_dir, go_prefix_enrich_xls)
                    copypdf(dest_dir, go_prefix_enrich_pdf)
                    go_prefix_enrich_xls,go_prefix_enrich_pdf = \
                        getRelativeDir([go_prefix_enrich_xls,go_prefix_enrich_pdf],
                                report_sub_dir)
                    go_prefix_enrich_xls_name,go_prefix_enrich_pdf_name = \
                        getFileName([go_prefix_enrich_xls,go_prefix_enrich_pdf])

                    print '''
(ref:enrich-{tag}-{count}-{enrich_cnt}) Top 30 gene ontology or KEGG enrichment terms for genes regulated by {type} in {samp}. 缺少的图表示对应样品中没有富集的条目。 Full lists of GO or KEGG enrichment terms can be downloaded ([{go_prefix_enrich_xls_name}]({go_prefix_enrich_xls})) ([{go_prefix_enrich_pdf_name}]({go_prefix_enrich_pdf})).             

```{{r enrich-{tag}-{count}-{enrich_cnt}, out.width="90%", fig.cap="(ref:enrich-{tag}-{count}-{enrich_cnt})"}}
knitr::include_graphics(c("{go_prefix_enrich_png}"))
```

'''.format(tag=tag2, count=count, enrich_cnt=enrich_cnt, type=type, samp=samp, 
        go_prefix_enrich_xls_name=go_prefix_enrich_xls_name, 
        go_prefix_enrich_pdf_name=go_prefix_enrich_pdf_name, 
        go_prefix_enrich_xls=go_prefix_enrich_xls, 
        go_prefix_enrich_pdf=go_prefix_enrich_pdf, 
        go_prefix_enrich_png=go_prefix_enrich_pdf.replace('pdf', 'png'))
#-------------------------------
#-------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    go_dir = options.go_dir
    sample_nameL = re.split(r'[, ]*', options.samp_name.strip())
    typeL = re.split(r'[, ]*', options.type.strip())
    gene_attrL = re.split(r'[, ]*', options.gene_attr.strip())
    enrich_type = options.enrich_type
    enrich_typeL = re.split(r'[, ]*', enrich_type.strip())
    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    verbose = options.verbose
    tag = options.tag
    global debug
    debug = options.debug
    #-----------------------------------
    aDict = {}
    generateDoc(report_dir, report_sub_dir, tag, go_dir, sample_nameL, typeL, gene_attrL, enrich_typeL)

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


