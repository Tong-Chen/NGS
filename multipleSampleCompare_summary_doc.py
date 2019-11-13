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
    This is designed to summarize output of multipleSampleCompare_summary.py, clusterProfileGo.sh, clusterprofileKEGG.sh. Mainly their **summary.xls** file will be parsed.


compare_summary:

[
    {
        "DE_count": {
            "file": "DE_genes/Pri_M2.all.DE.count.xls", 
            "plot": "DE_genes/Pri_M2.all.DE.count.xls.heatmapS.pdf"
        }, 
        "DE_parameters": {
            "fdr": 0.01, 
            "log2fc": 1.0
        }, 
        "DE_profile": {
            "file": "DE_genes/Pri_M2.all.DE.norm.kmeans.xls", 
            "plot": "DE_genes/Pri_M2.all.DE.norm.kmeans.sort.heatmapS.pdf"
        }, 
        "pca": {
            "plot": "Pri_M2.gene.xls.pca.scale.pdf"
        }, 
        "pearson": {
            "file": "Pri_M2.gene.xls.pearson.xls", 
            "plot": "Pri_M2.gene.xls.pearson.xls.pheatmap.pdf"
        }
    }, 
    {
        "Primodial\tM2": {
            "all": "DE_genes/Pri_M2.Primodial._vs_.M2.anno.xls", 
            "higherThan": "DE_genes/Pri_M2.Primodial._higherThan_.M2.anno.xls", 
            "higherThan_Cnt": 14, 
            "lowerThan": "DE_genes/Pri_M2.Primodial._lowerThan_.M2.anno.xls", 
            "lowerThan_Cnt": 2434
        }
    }
]

go-summary


{
    "Primodial\tM2": {
        "higherThan_BP": {
            "file": "DE_genes/Pri_M2.all.DE.entrez.Primodial._higherThan_.M2.BP_GO.xls", 
            "pdf": "DE_genes/Pri_M2.all.DE.entrez.Primodial._higherThan_.M2.BP_GO.scatterplot.dv.pdf"
        }, 
        "higherThan_CC": {
            "file": "DE_genes/Pri_M2.all.DE.entrez.Primodial._higherThan_.M2.CC_GO.xls", 
            "pdf": "DE_genes/Pri_M2.all.DE.entrez.Primodial._higherThan_.M2.CC_GO.scatterplot.dv.pdf"
        }, 
        "higherThan_MF": {
            "file": "DE_genes/Pri_M2.all.DE.entrez.Primodial._higherThan_.M2.MF_GO.xls", 
            "pdf": "DE_genes/Pri_M2.all.DE.entrez.Primodial._higherThan_.M2.MF_GO.scatterplot.dv.pdf"
        }, 
        "lowerThan_BP": {
            "file": "DE_genes/Pri_M2.all.DE.entrez.Primodial._lowerThan_.M2.BP_GO.xls", 
            "pdf": "DE_genes/Pri_M2.all.DE.entrez.Primodial._lowerThan_.M2.BP_GO.scatterplot.dv.pdf"
        }, 
        "lowerThan_CC": {
            "file": "DE_genes/Pri_M2.all.DE.entrez.Primodial._lowerThan_.M2.CC_GO.xls", 
            "pdf": "DE_genes/Pri_M2.all.DE.entrez.Primodial._lowerThan_.M2.CC_GO.scatterplot.dv.pdf"
        }, 
        "lowerThan_MF": {
            "file": "DE_genes/Pri_M2.all.DE.entrez.Primodial._lowerThan_.M2.MF_GO.xls", 
            "pdf": "DE_genes/Pri_M2.all.DE.entrez.Primodial._lowerThan_.M2.MF_GO.scatterplot.dv.pdf"
        }
    }
}

kegg_summary:

{
    "Primodial\tM2": {
        "higherThan_KEGG": {
            "file": "DE_genes/Pri_M2.all.DE.entrez.Primodial._higherThan_.M2.KEGG.xls", 
            "pdf": "DE_genes/Pri_M2.all.DE.entrez.Primodial._higherThan_.M2.KEGG.scatterplot.dv.pdf"
        }, 
        "lowerThan_KEGG": {
            "file": "DE_genes/Pri_M2.all.DE.entrez.Primodial._lowerThan_.M2.KEGG.xls", 
            "pdf": "DE_genes/Pri_M2.all.DE.entrez.Primodial._lowerThan_.M2.KEGG.scatterplot.dv.pdf"
        }
    }
}

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
import re
from tools import *
from json import load as json_load

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
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-c", "--compare-summary", dest="compare_summary",
        help="File name for summary.xls output by multipleSampleCompare_summary.py.")
    parser.add_option("-g", "--go-summary", dest="go_summary",
        help="File name for summary.xls output by clusterProfileGO.sh.")
    parser.add_option("-k", "--kegg-summary", dest="kegg_summary",
        help="File name for summary.xls output by clusterProfileKEGG.sh.")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='4_DE_gene_profile', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '4_DE_gene_profile'.")
    #parser.add_option("-d", "--doc-only", dest="doc_only",
    #    default=False, action="store_true", help="Specify to only generate doc.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.compare_summary != None, "A string needed for -c"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    compare_summary = options.compare_summary
    go_summary = options.go_summary
    kegg_summary = options.kegg_summary
    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)

    pic_label_base = grenerateLabel()

    print "# 样品基因表达差异比较 {#sample-compare}\n"
        
    sumL = json_load(open(compare_summary))

    if go_summary:
        goD = json_load(open(go_summary))
    else:
        goD = ''

    if kegg_summary:
        keggD = json_load(open(kegg_summary))
    else:
        keggD = ''
    
    # The first element if sumL containg 4 dicts with key as
    # DE_count, DE_profile, pca, pearson.
    # The value of these 4 keys is also a dict containing one or
    # two keys, 'file' and 'plot' which points to files.
    
    globalD = sumL[0]

    #pearson
    xls = globalD['pearson']['file']
    pdf = globalD['pearson']['plot']
    copy(dest_dir, xls)
    xls = getRelativeDir(xls, report_sub_dir)
    copypdf(dest_dir, pdf)
    pdf = getRelativeDir(pdf, report_sub_dir)
    png = pdf.replace('pdf', 'png')
    
    print """
**根据全基因表达量从从整体上对样品的相关性进行评估和样品分类，
确定相似性较高或较低的样品群 (Figure \@ref(fig:{pic_label_base}-pearson) and \@ref(fig:{pic_label_base}-pca))。**

(ref:{pic_label_base}-pearson) 样品间Pearson相关性分析。图中每个格子代表纵轴和横轴两个样品总体表达值的的Pearson correlation coefficient (PCC). 格子的颜色表示了PCC值的大小。[pdf]({pdf}) [XLS]({xls})

```{{r {pic_label_base}-pearson, fig.cap="(ref:{pic_label_base}-pearson)"}}
knitr::include_graphics("{png}")
```

""".format(pic_label_base=pic_label_base, png=png, pdf=pdf, xls=xls)

    #pca
    pdf = globalD['pca']['plot']
    copypdf(dest_dir, pdf)
    pdf = getRelativeDir(pdf, report_sub_dir)
    png = pdf.replace('pdf', 'png')
    
    print """

(ref:{pic_label_base}-pca) 样品主成分聚类 (PCA)分析。图中每个点代表一个样品，聚拢在一起的点具有相似的生物属性；相距越远的点代表样品的相似度越差。横轴表示第一主成分，对样品的区分度大于纵轴代表的第二主成分。[pdf]({pdf}) [XLS]({xls})

```{{r {pic_label_base}-pca, fig.cap="(ref:{pic_label_base}-pca)"}}
knitr::include_graphics("{png}")
```

""".format(pic_label_base=pic_label_base, png=png, pdf=pdf, xls=xls)
    
    #DE_parameters
    fdr = globalD['DE_parameters']['fdr']
    log2fc = globalD['DE_parameters']['log2fc']
    #DE_count
    xls = globalD['DE_count']['file']
    pdf = globalD['DE_count']['plot']
    copy(dest_dir, xls)
    xls = getRelativeDir(xls, report_sub_dir)
    copypdf(dest_dir, pdf)
    pdf = getRelativeDir(pdf, report_sub_dir)
    png = pdf.replace('pdf', 'png')

    print """
**样品间差异基因的鉴定和表达丰度展示 (FDR (false discovery rate) = {fdr}, log2FC (log2 fold change) = {log2fc})。**

样品间差异基因数目统计，每个方格代表横轴样品相比于纵轴样品丰度显著高的基因的数目 (Figure \@ref(fig:{pic_label_base}-de-count))。

(ref:{pic_label_base}-de-count) 样品间差异基因数目统计，每个方格代表横轴样品相比于纵轴样品丰度显著高的基因的数目。[PDF]({pdf}) [XLS]({xls})

```{{r {pic_label_base}-de-count, fig.cap="(ref:{pic_label_base}-de-count)"}}
knitr::include_graphics("{png}")
```
""".format(pic_label_base=pic_label_base, png=png, pdf=pdf, 
xls=xls, fdr=fdr, log2fc=log2fc)

    #DE_profile
    xls = globalD['DE_profile']['file']
    pdf = globalD['DE_profile']['plot']
    copy(dest_dir, xls)
    xls = getRelativeDir(xls, report_sub_dir)
    copypdf(dest_dir, pdf)
    pdf = getRelativeDir(pdf, report_sub_dir)
    png = pdf.replace('pdf', 'png')

    print """

样品间差异基因丰度图谱。每一行为一个基因，每一列为一个样品，每个小方格代表对应基因在对应样品的丰度，其值的大小根据颜色展示 (Figure \@ref(fig:{pic_label_base}-de-profile))。

(ref:{pic_label_base}-de-profile) 样品间差异基因丰度图谱。每一行为一个基因，每一列为一个样品，每个小方格代表对应基因在对应样品的丰度，其值的大小根据颜色展示。[PDF]({pdf}) [XLS]({xls})

```{{r {pic_label_base}-de-profile, fig.cap="(ref:{pic_label_base}-de-profile)"}}
knitr::include_graphics("{png}")
```
""".format(pic_label_base=pic_label_base, png=png, pdf=pdf, 
xls=xls, fdr=fdr, log2fc=log2fc)

    # DE compare -- each sample pair

    # The second element of sumL is a dict containing sample-compare 
    # information in sub-dict.
    compD = sumL[1]

    count = 0
    for comp_pair, comp_subD in compD.items():
        count += 1
        condA, condB = comp_pair.split('\t')
        all = comp_subD['all']
        higherThan = comp_subD['higherThan']
        higherThan_Cnt = comp_subD['higherThan_Cnt']
        lowerThan = comp_subD['lowerThan']
        lowerThan_Cnt = comp_subD['lowerThan_Cnt']
        copy(dest_dir, all, higherThan, lowerThan)
        all, higherThan, lowerThan = \
            getRelativeDir([all, higherThan, lowerThan], report_sub_dir)
        print """
#### Sample compare ({condA} vs {condB}) (bin_size={bin})

{condA}和{condB}检测到的基因Reads丰度统计和基因注释 ([点击下载)](all)。如果对应基因落在基因启动子基因 (转录起始位点上游1 kb和下游500 bp基因), 则以基因的注释作为此基因注释。一个基因可能对应多个基因。

{condA}相比于{condB}显著高的基因有{higherThan_Cnt}个，注释见[表格]({higherThan}).

{condB}相比于{condA}显著高的基因有{lowerThan_Cnt}个，注释见[表格]({lowerThan}).

""".format(condA=condA, condB=condB, bin=bin, all=all, higherThan_Cnt=higherThan_Cnt, lowerThan_Cnt=lowerThan_Cnt, higherThan=higherThan, lowerThan=lowerThan)

        if keggD:
            kegg_annoD = keggD[comp_pair]
            higherThan_kegg_file = kegg_annoD['higherThan_KEGG']['file']
            higherThan_kegg_pdf = kegg_annoD['higherThan_KEGG']['pdf']
            higherThan_kegg_png = higherThan_kegg_pdf.replace('pdf', 'png')
            copy(dest_dir, higherThan_kegg_file)
            copypdf(dest_dir, higherThan_kegg_pdf)
            higherThan_kegg_file, higherThan_kegg_pdf, higherThan_kegg_png = \
                getRelativeDir(\
                [higherThan_kegg_file, higherThan_kegg_pdf, higherThan_kegg_png], 
                report_sub_dir)
            lowerThan_kegg_file = kegg_annoD['lowerThan_KEGG']['file']
            lowerThan_kegg_pdf = kegg_annoD['lowerThan_KEGG']['pdf']
            lowerThan_kegg_png = lowerThan_kegg_pdf.replace('pdf', 'png')
            copy(dest_dir, lowerThan_kegg_file)
            copypdf(dest_dir, lowerThan_kegg_pdf)
            lowerThan_kegg_file, lowerThan_kegg_pdf, lowerThan_kegg_png = \
                getRelativeDir(\
                [lowerThan_kegg_file, lowerThan_kegg_pdf, lowerThan_kegg_png], 
                report_sub_dir)

            print '''
(ref:{pic_label_base}-{count}-kegg) Top 30 KEGG enrichment terms for genes up-regulated in {condA} (left) and {condB} (right). 缺少的图表示对应样品中没有富集的条目。 Full lists of KEGG enrichment terms can be downloaded for **{condA}** ([Enrich Table]({higherThan_kegg_file})) ([ PDF pic]({higherThan_kegg_pdf})) and **{condB}** ([Enrich Table]({lowerThan_kegg_file})) ([ PDF pic]({lowerThan_kegg_pdf})).             

```{{r {pic_label_base}-{count}-kegg, out.width="50%", fig.cap="(ref:{pic_label_base}-{count}-kegg)"}}
knitr::include_graphics(c("{higherThan_kegg_png}", "{lowerThan_kegg_png}"))
```

'''.format(pic_label_base=pic_label_base, count=count, condA=condA, condB=condB, 
        higherThan_kegg_file=higherThan_kegg_file, 
        higherThan_kegg_pdf=higherThan_kegg_pdf, 
        higherThan_kegg_png=higherThan_kegg_png, 
        lowerThan_kegg_file=lowerThan_kegg_file, 
        lowerThan_kegg_pdf=lowerThan_kegg_pdf, 
        lowerThan_kegg_png=lowerThan_kegg_png)

        if goD:
            go_annoD = goD[comp_pair]
            goL = ['BP', 'CC', 'MF']
            for go in goL:
                higherThan_go_file = go_annoD['higherThan_'+go]['file']
                higherThan_go_pdf = go_annoD['higherThan_'+go]['pdf']
                higherThan_go_png = higherThan_go_pdf.replace('pdf', 'png')
                copy(dest_dir, higherThan_go_file)
                copypdf(dest_dir, higherThan_go_pdf)
                higherThan_go_file, higherThan_go_pdf, higherThan_go_png = getRelativeDir([higherThan_go_file, higherThan_go_pdf, higherThan_go_png], report_sub_dir)
                lowerThan_go_file = go_annoD['lowerThan_'+go]['file']
                lowerThan_go_pdf = go_annoD['lowerThan_'+go]['pdf']
                lowerThan_go_png = lowerThan_go_pdf.replace('pdf', 'png')
                copy(dest_dir, lowerThan_go_file)
                copypdf(dest_dir, lowerThan_go_pdf)
                lowerThan_go_file, lowerThan_go_pdf, lowerThan_go_png = getRelativeDir([lowerThan_go_file, lowerThan_go_pdf, lowerThan_go_png], report_sub_dir)

                print '''
(ref:{pic_label_base}-{count}-{go}) Top 30 Gene ontology ({go}) enrichment terms for genes up-regulated in {condA} (left) and {condB} (right). 缺少的图表示对应样品中没有富集的条目。 Full lists of Gene ontology ({go}) enrichment terms can be downloaded for **{condA}** ([Enrich Table]({higherThan_go_file})) ([ PDF pic]({higherThan_go_pdf})) and **{condB}** ([Enrich Table]({lowerThan_go_file})) ([ PDF pic]({lowerThan_go_pdf})).             

```{{r {pic_label_base}-{count}-{go}, out.width="50%", fig.cap="(ref:{pic_label_base}-{count}-{go})"}}
knitr::include_graphics(c("{higherThan_go_png}", "{lowerThan_go_png}"))
```

'''.format(pic_label_base=pic_label_base, go=go, count=count, 
        condA=condA, condB=condB, 
        higherThan_go_file=higherThan_go_file, 
        higherThan_go_pdf=higherThan_go_pdf, 
        higherThan_go_png=higherThan_go_png, 
        lowerThan_go_file=lowerThan_go_file, 
        lowerThan_go_pdf=lowerThan_go_pdf, 
        lowerThan_go_png=lowerThan_go_png)


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


