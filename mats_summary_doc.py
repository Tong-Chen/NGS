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
    This is designed to get the result of MATS.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from tools import *
import re

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
    parser.add_option("-c", "--compare-pair", dest="filein",
        metavar="COMPARE_PAIR", 
        help="compare_pair file given to MATS.py")
    parser.add_option("-p", "--prefix", dest="prefix",
        metavar="PREFIX", 
        help="The prefix for input files.")
    parser.add_option("-q", "--fdr", dest="fdr",
        default="0.05", help="FDR for selecting significant AS")
    parser.add_option("-g", "--go-prefix", dest="go_prefix",
        help="String like <MATS/go/$(prefix).mats.allType.DE.entrez>.")
    parser.add_option("-e", "--enrichmentType", dest="enrich_type",
        default="BP_GO, MF_GO, CC_GO, KEGG", 
        help="Default <BP_GO, MF_GO, CC_GO, KEGG>. And only these 4 types of enrichment are supported. These strings will be used to get filenames.")
    #parser.add_option("-o", "--output-file-prefix", dest="op",
    #    help="Specify output file prefix")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='4_alternative_splicing', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '4_alternative_splicing'.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    compare_pair = options.filein
    prefix = options.prefix
    fdr    = options.fdr
    verbose = options.verbose
    go_prefix = options.go_prefix
    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    global debug
    debug = options.debug
    enrich_type = options.enrich_type
    enrich_typeL = re.split(r'[, ]*', enrich_type.strip())
    #-----------------------------------
    typeL = ['SE', 'MXE', 'A5SS', 'A3SS', 'RI']
    nameD = {'SE':'Exon skipping', 'MXE':'Mutually exclusive exon', 
            'RI':'Intron retention', 
            'A5SS':"Alternative 5' splice site", 
            'A3SS':"Alternative 3' splice site"}

    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)
    result_dir = 'MATS/'

    print "# 选择性剪接分析\n"
    prefix = result_dir + prefix
    as_all_file = prefix + ".AS.jcCount"
    as_all = as_all_file + ".stackBars.pdf"
    as_sig_all_file = prefix + ".AS.sigJcCount"
    as_sig_all = as_sig_all_file + ".stackBars.pdf"
    as_sig_all_up_dw_file = prefix + ".AS.sigJcCount.up_dw"
    as_sig_all_up_dw = as_sig_all_up_dw_file + ".stackBars.pdf"
    copy(dest_dir, as_all_file, as_sig_all_file, as_sig_all_up_dw_file)
    as_all_file, as_sig_all_file, as_sig_all_up_dw_file = getRelativeDir([as_all_file,  as_sig_all_file,  as_sig_all_up_dw_file], report_sub_dir)

    copypdf(dest_dir, as_all, as_sig_all, as_sig_all_up_dw)
    as_all, as_sig_all, as_sig_all_up_dw = getRelativeDir([as_all, as_sig_all, as_sig_all_up_dw],report_sub_dir)

    copy(report_dir+'/images/', "/MPATHB/self/resource/sample/splicing.png", "/MPATHB/self/resource/sample/MATS_explanation.png")
    print '''
## 选择性剪接总览
'''
    knitr_read_txt(report_dir, "DE_AS_summary")

    print '''
可变剪接事件主要有5种，如 Figure \@ref(fig:as-event-example)所示。Figure \@ref(fig:as-event-example)右图展示了基本的可变剪接分析思路和算法演示。

```{r as-event-example, fig.cap="Illustration of five types of alternative splicing events.[ref](http://rnaseq-mats.sourceforge.net/)"}
knitr::include_graphics("images/splicing.png")
```

样品之间两两比较鉴定出的可变剪接事件的分布如 Figure \@ref(fig:all-as-event), 样品间显著差异(fdr<%s)的可变剪接事件分布如 Figure \@ref(fig:de-as-event)和\@ref(fig:de-as-event-2)所示。

''' % fdr

    print '''
(ref:all-as-event) Statistics of all AS events between each two samples. [XLS]({xls}) [PDF]({pdf})

```{{r all-as-event, fig.cap="(ref:all-as-event)" }}
knitr::include_graphics("{png}")
```
'''.format(xls=as_all_file, pdf=as_all, png=as_all.replace('pdf', 'png'))

    print '''
(ref:de-as-event) Statistics of all significantly differential AS events between each two samples. [XLS]({xls}) [PDF]({pdf})

```{{r de-as-event, fig.cap="(ref:de-as-event)"}}
knitr::include_graphics("{png}")
```
'''.format(xls=as_sig_all_file, pdf=as_sig_all, png=as_sig_all.replace('pdf', 'png'))

    print '''
(ref:de-as-event-2) Statistics of all significantly differential AS events between each two samples. [XLS]({xls}) [PDF]({pdf})

```{{r de-as-event-2, fig.cap="(ref:de-as-event-2)" }}
knitr::include_graphics("{png}")
```
'''.format(xls=as_sig_all_up_dw_file, pdf=as_sig_all_up_dw, png=as_sig_all_up_dw.replace('pdf', 'png'))

    count = 0
    print "\n## 样品间选择性剪接事件比较\n"

    print "**结果表格中表头部分的解释 Table \@ref(tab:table-header-mats-output) and Figure \@ref(fig:mats-as-event-header-explanation).**"
    
    print '''

Table: (\#tab:table-header-mats-output) 选择性剪接分析和注释表格中表头解释。

```{r table-header-mats-output}
table_header_mats_output <- "表头;解释
ID;可变剪接分析编号，标示使用，无特殊含义
GeneID;注释的基因的名字
geneSymbol;基因名字
chr;染色体名字
strand;基因处在染色体的正链还是负链
riExonStart_0base;保留内含子后的新外显子在基因组的起始位置 (0-based, included)
riExonEnd;保留内含子后的新外显子在基因组的终止位置 (0-based, not included)
upstreamES;发生可变剪接的外显子上游外显子起始位置 (0-based, included)
upstreamEE;发生可变剪接的外显子上游外显子终止位置 (0-based, not included)
downstreamES;发生可变剪接的外显子下游外显子起始位置 (0-based, included)
downstreamEE;发生可变剪接的外显子下游外显子终止位置 (0-based, not included)
longExonStart_0base;A3SS中下游外显子最左侧剪接位点或A5SS中上游外显子起始
longExonEnd;A3SS中下游外显子终止位置或A5SS中上游外显子最右侧剪接位点
shortES;A3SS中下游外显子次左侧剪接位点或A5SS中上游外显子起始
shortEE;A3SS中下游外显子终止位置或A5SS中上游外显子次右侧剪接位点
flankingES;A3SS中上游外显子起始位置或A5SS中下游外显子起始位置
flankingEE;A3SS中上游外显子终止位置或A5SS中下游外显子终止位置
exonStart;SE中保留的外显子的起始位置
exonEnd;SE中保留的外显子的终止位置
1stExonStart_0base;MXE中第一个使用的外显子的起始位置
1stExonEnd;MXE中第一个使用的外显子的终止位置
2ndExonStart_0base;MXE中第二个使用的外显子的起始位置
2ndExonEnd;MXE中第二个使用的外显子的终止位置
IJC_T40;T40样品中包含可变剪接位置处Junction reads数，不同的重复用逗号分开
SJC_T40;T40样品中去除可变剪接位置处Junction reads数，不同的重复用逗号分开
IJC_T45;T45样品中包含可变剪接位置处Junction reads数，不同的重复用逗号分开
SJC_T45;T45样品中去除可变剪接位置处Junction reads数，不同的重复用逗号分开
IncFormLen;包含可变剪接部分的转录本的长度
SkipFormLen;不包含可变剪接部分的转录本长度
Pvalue;差异检验显著性度量
FDR;多重假设检验校正的假阳性率
IncLevel1;样品1（T40)中包含可变剪接部分的转录本的表达，不同的重复用逗号分开
IncLevel2;样品2（T45)中包含可变剪接部分的转录本的表达，不同的重复用逗号分开
IncLevelDifference;样品1(T40)的转录本表达量平均值减去样品2的转录本表达平均值(T45)
Remaining columns;Function annotation"

table_header_mats_output_mat <- read.table(text=table_header_mats_output, header=T, sep=";")

knitr::kable(table_header_mats_output_mat, boktabs=T, format="markdown")
```

```{r mats-as-event-header-explanation, fig.cap="Illustration of 5 types AS events and explanation of MATS result. ES: exon start; EE: exon end. "}
knitr::include_graphics("images/MATS_explanation.png")
```

'''

    for line in open(compare_pair):
        samp1, samp2 = line.split()
        count += 1
        print "### 样品 (%s vs %s) 选择性剪接事件分析\n" % (samp1,samp2)
        samp = samp1+'_vs_'+samp2

        as_eventL = ["AS event;Acronym;File"]
        for type in typeL:
            type_f = '.'.join([samp, type, 'MATS.JunctionCountOnly.txt.anno.xls'])
            copy(dest_dir, result_dir+type_f)
            tmp = [nameD[type], type, "[{}]({})".format(type_f, report_sub_dir+'/'+type_f)]
            as_eventL.append(';'.join(tmp))

        print '''
Table: (\#tab:{label}) List of annotated total AS events.          

```{{r {label} }}
{label} <- "{table_context}"
{label}_mat <- read.table(text={label}, sep=";", header=T, quote="")
knitr::kable({label}_mat, booktabs=T, format="markdown")
```
'''.format(label="asEventTotal"+str(count), table_context='\n'.join(as_eventL))

        #print "* Annotation of %s up-regulated AS events\n" % samp1
        as_eventL = ["AS event;Acronym;Type;File"]
        samp1_up = samp + '.' + samp1 + '.UP.' + fdr
        for type in typeL:
            type_f = '.'.join([samp1_up, type, 'MATS.JunctionCountOnly.txt.anno.xls'])
            copy(dest_dir, result_dir+type_f)
            tmp = [nameD[type], type, samp1+'.UP', "[{}]({})".format(type_f, report_sub_dir+'/'+type_f)]
            as_eventL.append(';'.join(tmp))

        #print "* Annotation of %s up-regulated AS events\n" % samp2
        samp2_up = samp + '.' + samp2 + '.UP.' + fdr
        for type in typeL:
            type_f = '.'.join([samp2_up, type, 'MATS.JunctionCountOnly.txt.anno.xls'])
            copy(dest_dir, result_dir+type_f)
            tmp = [nameD[type], type, samp2+'.UP', "[{}]({})".format(type_f, report_sub_dir+'/'+type_f)]
            as_eventL.append(';'.join(tmp))

        print '''
Table: (\#tab:{label}) List of significantly differentiated AS events between {samp1} and {samp2}.          

```{{r {label} }}
{label} <- "{table_context}"
{label}_mat <- read.table(text={label}, sep=";", header=T, quote="")
knitr::kable({label}_mat, booktabs=T, format="markdown")
```
'''.format(label="asEventTotalSig"+str(count), table_context='\n'.join(as_eventL), 
    samp1=samp1, samp2=samp2 )
        
        print "\n#### 差异选择性剪接基因的功能注释\n\n"
        go_prefix_up = go_prefix+'.'+samp+'.'+samp1+'.UP.'
        go_prefix_dw = go_prefix+'.'+samp+'.'+samp2+'.UP.'
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
(ref:AS-enrich-{count}-{enrich_cnt}) Top 30 gene ontology or KEGG enrichment terms for alternative spliced genes up-regulated in {samp1} (left) and {samp2} (right). 缺少的图表示对应样品中没有富集的条目。 Full lists of GO or KEGG enrichment terms can be downloaded for **{samp1}** ([Enrich Table]({rgo_prefix_up_enrich_xls})) ([ PDF pic]({rgo_prefix_up_enrich_pdf})) and **{samp2}** ([Enrich Table]({rgo_prefix_dw_enrich_xls})) ([ PDF pic]({rgo_prefix_dw_enrich_pdf})).             

```{{r AS-enrich-{count}-{enrich_cnt}, out.width="50%", fig.cap="(ref:AS-enrich-{count}-{enrich_cnt})"}}
knitr::include_graphics(c("{rgo_prefix_up_enrich_png}", "{rgo_prefix_dw_enrich_png}"))
```

'''.format(count=count, enrich_cnt=enrich_cnt, samp1=samp1, samp2=samp2, rgo_prefix_up_enrich_xls=rgo_prefix_up_enrich_xls, rgo_prefix_up_enrich_pdf=rgo_prefix_up_enrich_pdf, rgo_prefix_dw_enrich_xls=rgo_prefix_dw_enrich_xls, rgo_prefix_dw_enrich_pdf=rgo_prefix_dw_enrich_pdf, rgo_prefix_up_enrich_png=rgo_prefix_up_enrich_png, rgo_prefix_dw_enrich_png=rgo_prefix_dw_enrich_png)
#-------------------------------

        #print "* Plots for AS events\n"
        print "\n#### 差异选择性剪接基因图示\n\n"

        as_eventL = ["AS event;Acronym;Type;Plot"]
        for type in typeL:
            type_f = '.'.join([samp1_up, type, 'MATS.JunctionCountOnly.sashimi'])
            type_f_zip = type_f+'.zip'       
            os.system("zip %s/%s -r %s/%s/Sashimi_plot >/dev/null 2>&1" \
                % (dest_dir, type_f_zip, result_dir, type_f))
            tmp = [nameD[type], type, samp1+'.UP', "[{}]({})".format(type_f, report_sub_dir+'/'+type_f_zip)]
            as_eventL.append(';'.join(tmp))

        for type in typeL:
            type_f = '.'.join([samp2_up, type, 'MATS.JunctionCountOnly.sashimi'])
            type_f_zip = type_f+'.zip'       
            os.system("zip %s/%s -r %s/%s/Sashimi_plot >/dev/null 2>&1" \
                % (dest_dir, type_f_zip, result_dir, type_f))
            tmp = [nameD[type], type, samp2+'.UP', "[{}]({})".format(type_f, report_sub_dir+'/'+type_f_zip)]
            as_eventL.append(';'.join(tmp))

        print '''
Table: (\#tab:{label}) Sashimi-plot of significantly differentiated AS events between {samp1} and {samp2}.          

```{{r {label} }}
{label} <- "{table_context}"
{label}_mat <- read.table(text={label}, sep=";", header=T, quote="")
knitr::kable({label}_mat, booktabs=T, format="markdown")
```
'''.format(label="asEventTotalSigPlot"+str(count), table_context='\n'.join(as_eventL), 
    samp1=samp1, samp2=samp2 )

        print "\n**Examples for each type of AS events**\n"
        for type in typeL:
            type_f1 = '.'.join([result_dir+samp1_up, type, \
                'MATS.JunctionCountOnly.sashimi/Sashimi_plot'])
            type_f2 = '.'.join([result_dir+samp2_up, type, \
                'MATS.JunctionCountOnly.sashimi/Sashimi_plot'])
            #type_f2 = '.'.join([samp2_up, type, 'MATS.JunctionCountOnly.sashimi'])
            #pdf1 = pdf2 = png1 = png2 = 'unexist'
            if not os.path.exists(type_f1):
                os.system("mkdir -p "+type_f1)
                os.system("/bin/cp -u /MPATHB/self/resource/sample/no_result.pdf "+type_f1)
            type_f1_l = os.listdir(type_f1)
            tmp_f_1 = type_f1 + '/' + type_f1_l[0]
            copypdf(dest_dir, tmp_f_1)
            pdf1=report_sub_dir+'/'+type_f1_l[0]
            png1 = pdf1.replace('pdf', 'png')
            if not os.path.exists(type_f2):
                os.system("mkdir -p "+type_f2)
                os.system("/bin/cp -u /MPATHB/self/resource/sample/no_result.pdf "+type_f2)
            type_f2_l = os.listdir(type_f2)
            tmp_f_2 = type_f2 + '/' + type_f2_l[0]
            copypdf(dest_dir, tmp_f_2)
            pdf2=report_sub_dir+'/'+type_f2_l[0]
            png2 = pdf2.replace('pdf', 'png')

            print '''
(ref:{label}) Sashimiplot for **{type}** up-regulated in {samp1} (left, [PDF]({pdf1})) and {samp2} (right, [PDF]({pdf2})) (未显示的图形表示对应类型的剪接在对应样品没有差异).

```{{r {label}, out.width="49%", fig.cap="(ref:{label})" }}
knitr::include_graphics(c("{png1}", "{png2}"))
```
'''.format(label=type+'sig'+str(count), type=type, samp1=samp1, samp2=samp2, 
    pdf1=pdf1, pdf2=pdf2, png1=png1, png2=png2)
        
    #-----------------------------------
       
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


