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
        metavar="FILEIN", help="The all.DE file generated using DESeq.sh. ")
    parser.add_option("-e", "--enrichmentType", dest="enrich_type",
        default="BP_GO, MF_GO, CC_GO, KEGG", 
        help="Default <BP_GO, MF_GO, CC_GO, KEGG>. And only these 4 types of enrichment are supported. These strings will be used to get filenames.")
    parser.add_option("-g", "--go-prefix", dest="go_prefix",
        help="String like <go/NK_trans.rc.xls.DESeq2.all.DE.entrez>.")
    parser.add_option("-F", "--log2-FC", dest="log2FC",
        help="log2_FC for screening DE genes")
    parser.add_option("-q", "--FDR", dest="FDR",
        help="FDR for screening DE genes")
    parser.add_option("-t", "--report-tag", dest="report_tag",
        default='de', help="A unique string, only alphabets allowed. Default 'de'.")
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

def generateDoc(report_dir, report_sub_dir, report_tag, file, go_prefix, enrich_typeL, curation_label, log2FC, FDR):
    
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)

    print "\n## 基因表达图谱和差异基因 {#expr-profile-de-genes}\n"
    
    
    print "\n### 整体表达聚类"

    knitr_read_txt(report_dir, "DE_summary_cluster_pca")

    prefix = file.replace('.all.DE', '')
    '''
    ## OLD
    #typeL = ['A._vs_.C_up', 'A._vs_.C_dw']
    typeL = ['M1._higherThan_.Y1', 'M1._lowerThan_.Y1']
    '''
    typeL = set([line.split()[1] for line in open(file)])
    '''
    typeD = {'A._vs_.C':['A._higherThan_.C', 'A._lowerThan_.C']}
    '''
    typeD = {}
    for cmp_type in typeL:
        firstSep = cmp_type.find('._')
        secndSep = cmp_type.find('_.')
        sampA = cmp_type[:firstSep]
        sampB = cmp_type[secndSep+2:]
        type = cmp_type[firstSep:secndSep+2]
        #name = cmp_type.replace('._higherThan_.', '____')
        #name = name.repalce('._lowerThan_.', '____')
        #nameL = name.split('___')
        cmp = '._vs_.'.join([sampA, sampB])
        #cmp, type = cmp_type.rsplit('_', 1)
        type = type.join([sampA, sampB])
        if cmp not in typeD:
            typeD[cmp] = [type]
        else:
            typeD[cmp].append(type)
    typeL = typeD.keys()        
    #--------------------------------------------        
    anno_all = prefix + '.normalized.anno.xls'
    pca = prefix + ".normalized.rlog.pca.pdf"
    pearson = prefix + ".normalized.rlog.pearson.pdf"
    pca_png = prefix + ".normalized.rlog.pca.png"
    pearson_png = prefix + ".normalized.rlog.pearson.png"

    copy(dest_dir, anno_all)
    copypdf(dest_dir, pca, pearson)
    
    rAnno, rPca, rPearson = getRelativeDir([anno_all, pca_png, pearson_png], 
            report_sub_dir)
    print '''
样品表达谱的相关性分析和主成分分析，展示的是样品之间基因表达的整体相似性和差异以及样品之间的聚类关系。

所有样品基因表达及注释文件下载：[基因表达和注释矩阵]({rAnno})

(ref:sample-pca-pearson) Pearson correlation and PCA analysis for all samples. 左图首先计算样品两两之间基因表达的Pearson相关系数，然后层级聚类分析样品之间的相似性，最后绘制Heatmap展示。每个方块的颜色代表横轴和纵轴对应样品的表达谱的相似度，值越高表明两个样品的基因表达谱越相近。右图是转换原始表达数据为两个正交的主成分，然后绘制样品在这两个主成分组成的空间内的分布关系，空间关系越近的样品相似度越高。

```{{r sample-pca-pearson, out.width="50%", fig.cap="(ref:sample-pca-pearson)"}}
knitr::include_graphics(c("{rPearson}", "{rPca}"))
```
'''.format(rAnno=rAnno, rPearson=rPearson, rPca=rPca)
    
    print "\n### 差异基因功能注释和富集分析\n"
    
    knitr_read_txt(report_dir, "DE_summary_gene_anno_enrichment")

    print "\n#### 差异基因总结及表达谱\n"
    
    DE_count = file + ".count.xls"
    copy(dest_dir, DE_count)
    DE_count_pdf = DE_count + '.heatmapS.pdf'
    copypdf(dest_dir, DE_count_pdf)
    DE_count_png = DE_count_pdf.replace('pdf', 'png')

    DE_count, DE_count_pdf, DE_count_png = \
        getRelativeDir([DE_count, DE_count_pdf, DE_count_png], report_sub_dir)
    label_tmp_146 = 'de-count-heat-146'

    print """
样品间差异基因数目统计，每个方格代表横轴样品相比于纵轴样品丰度显著高的基因的数目 (Figure \@ref(fig:{label}))。

(ref:{label}) 样品间差异基因数目统计，每个方格代表横轴样品相比于纵轴样品丰
度显著高的基因的数目。白色方格表示对应样品未做比较。[PDF]({pdf}) [XLS]({xls})

```{{r {label}, fig.cap="(ref:{label})"}}
knitr::include_graphics("{png}")
```
""".format(label=label_tmp_146+report_tag, pdf=DE_count_pdf, png=DE_count_png, xls=DE_count)

    DE_profile = file + ".norm.kmeans.xls"
    copy(dest_dir, DE_profile)
    DE_profile_pdf = file + '.norm.kmeans.sort.heatmapS.pdf'
    copypdf(dest_dir, DE_profile_pdf)
    DE_profile_png = DE_profile_pdf.replace('pdf', 'png')

    DE_profile, DE_profile_pdf, DE_profile_png = \
        getRelativeDir([DE_profile, DE_profile_pdf, DE_profile_png], report_sub_dir)
    label_tmp_146 = 'de-profile-heat-146'+report_tag

    print """
样品间差异基因表达图谱，每一行代表一个基因，每一列代表一组样品，每个小方格代表基因在对应样品的相对表达量，其值的大小用颜色表示 (Figure \@ref(fig:{label}))。

(ref:{label}) 样品间差异基因表达图谱，每一行代表一个基因，每一列代表一组样品，每个小方格代表基因在对应样品的相对表达量，其值的大小用颜色表示。[PDF]({pdf}) [XLS]({xls})

```{{r {label}, fig.cap="(ref:{label})"}}
knitr::include_graphics("{png}")
```
""".format(label=label_tmp_146, pdf=DE_profile_pdf, png=DE_profile_png, xls=DE_profile)

    count = 0
    '''
    typeD = {'A._vs_.C':['A._higherThan_.C', 'A._lowerThan_.C']}
    '''
    for type in typeL:
        count += 1
        volcano = prefix + '.'+type+'.Volcano.png'
        venn = file + '.'+type+'.vennDiagram.pdf'
        venn_png = venn.replace('pdf', 'png')
        copypng(dest_dir, volcano)
        copypdf(dest_dir, venn)
        samp1, samp2 = type.split("._vs_.")
        anno_up_de = prefix+'.'+'._higherThan_.'.join([samp1, samp2])+'.xls.anno.xls'
        anno_dw_de = prefix+'.'+'._lowerThan_.'.join([samp1, samp2])+'.xls.anno.xls'
        #anno_dw_de = prefix+'.'+type+'.results.DE_dw.anno.xls'
        copy(dest_dir, anno_up_de, anno_dw_de)

        ranno_up_de, ranno_dw_de, rvolcano, rvenn = \
            getRelativeDir([anno_up_de, anno_dw_de, volcano, venn_png], report_sub_dir)

        print "\n#### DE genes between {}\n".format(type.replace("._vs_.", " and "))
        print """

(ref:volcano-venn-{count}) Volcano plot and venn diagram showing DE genes (FDR<{FDR}, |log2FC|>={log2FC}). 左图火山图展示的是基因差异倍数和显著性FDR的关系，红点表示鉴定出的差异表达基因。log2FC大于0的部分为在样品**{samp1}**中高表达的基因；log2FC小于0的部分为在样品**{samp2}**中高表达的基因。右图Venn图展示的是上调和下调的基因数目。{samp1}._vs_.{samp2}.up表示{samp1}相比于{samp2}高表达的基因。

```{{r volcano-venn-{count}, out.width="50%", fig.cap="(ref:volcano-venn-{count})"}}
knitr::include_graphics(c("{volcano}", "{venn}"))
```
""".format(count=str(count)+report_tag, volcano=rvolcano, venn=rvenn, samp1=samp1, samp2=samp2, FDR=FDR, log2FC=log2FC)
        
        print """
* 样品{samp1}上调的基因及其注释: [Full table]({anno_up_de})
* 样品{samp2}上调的基因及其注释: [Full table]({anno_dw_de})
""".format(samp1=samp1, samp2=samp2, anno_up_de=ranno_up_de, anno_dw_de=ranno_dw_de)

        print """ 
(ref:anno-up-de-{count}) 样品**{samp1}**上调的基因及其注释(top 30, first 7 columns)。

```{{r anno-up-de-{count} }}
anno_up_de <- read.table("{anno_up_de}",header=T,row.names=1, nrows=31, quote="",check.names=F, sep="\\t", comment.char="")
knitr::kable(head(anno_up_de, 30)[, 1:7], booktabs=T, caption="(ref:anno-up-de-{count})")
```

(ref:anno-dw-de-{count}) 样品**{samp2}**上调的基因及其注释(top 30, first 7 columns)。

```{{r anno-dw-de-{count} }}
anno_dw_de <- read.table("{anno_dw_de}",header=T,row.names=1, nrows=31, quote="",check.names=F, sep="\\t", comment.char="")
knitr::kable(head(anno_dw_de, 30)[, 1:7], booktabs=T, caption="(ref:anno-dw-de-{count})")
```

""".format(count=str(count)+report_tag, samp1=samp1, samp2=samp2, anno_up_de=ranno_up_de,  anno_dw_de=ranno_dw_de) 

        print '''
Gene ontology 和 KEGG富集结果展示。GO富集分3个门类，生物进程 (BP)、分子功能 (MF)、细胞组分 (CC)。具体见下方对应的图和图例解释。蓝色字体为超链接，点击可下载。若有某一门类缺失，表示无富集或无对应结果。
        
        '''

        '''
        typeD = {'A._vs_.C':['A._higherThan_.C', 'A._lowerThan_.C']}
        '''
        go_prefix_up = go_prefix+'.'+'._higherThan_.'.join([samp1,samp2])+'.'
        go_prefix_dw = go_prefix+'.'+'._lowerThan_.'.join([samp1,samp2])+'.'
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
(ref:enrich-{count}-{enrich_cnt}) Top 30 gene ontology or KEGG enrichment terms for genes up-regulated in {samp1} (left) and {samp2} (right). Nodes color green to red represents enrichment signifiance. Node sizes represents number of DE genes in given terms. 点的颜色代表在对应通路的显著性富集程度，颜色越红代表富集越显著。点的大小代表对应通路中差异基因的数目，点越大，对应通路差异基因越多。缺少的图表示对应样品中没有富集的条目。 Full lists of GO or KEGG enrichment terms can be downloaded for **{samp1}** ([Enrich Table]({rgo_prefix_up_enrich_xls})) ([ PDF pic]({rgo_prefix_up_enrich_pdf})) and **{samp2}** ([Enrich Table]({rgo_prefix_dw_enrich_xls})) ([ PDF pic]({rgo_prefix_dw_enrich_pdf})).             

```{{r enrich-{count}-{enrich_cnt}, out.width="50%", fig.cap="(ref:enrich-{count}-{enrich_cnt})"}}
knitr::include_graphics(c("{rgo_prefix_up_enrich_png}", "{rgo_prefix_dw_enrich_png}"))
```

'''.format(count=str(count)+report_tag, enrich_cnt=enrich_cnt, samp1=samp1, samp2=samp2, rgo_prefix_up_enrich_xls=rgo_prefix_up_enrich_xls, rgo_prefix_up_enrich_pdf=rgo_prefix_up_enrich_pdf, rgo_prefix_dw_enrich_xls=rgo_prefix_dw_enrich_xls, rgo_prefix_dw_enrich_pdf=rgo_prefix_dw_enrich_pdf, rgo_prefix_up_enrich_png=rgo_prefix_up_enrich_png, rgo_prefix_dw_enrich_png=rgo_prefix_dw_enrich_png)
#-------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    enrich_type = options.enrich_type
    enrich_typeL = re.split(r'[, ]*', enrich_type.strip())
    log2_FC = options.log2FC
    FDR = options.FDR
    report_tag = options.report_tag
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
    generateDoc(report_dir, report_sub_dir, report_tag, file, go_prefix, enrich_typeL, curation_label, log2_FC, FDR)

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


