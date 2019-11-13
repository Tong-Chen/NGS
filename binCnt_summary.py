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
    This is designed to summarize output of meltRegionBinMatrix.py, multipleSampleCompare_summary.py. Mainly their **summary.xls** file will be parsed.
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
        default="binC0304YCnt.merge/SCX.binC0304YCnt.merge.all.DE.summary.xls", 
        help="File name for summary.xls output by multipleSampleCompare_summary.py. One special symbol <C0304Y> should be put at the bin-size position. For example, file <bin10000Cnt.merge/SCX.bin10000Cnt.merge.all.DE.summary.xls> should be <binC0304YCnt.merge/SCX.binC0304YCnt.merge.all.DE.summary.xls> (which is default). The program will replace <C0304Y> with given bin sizes.")
    parser.add_option("-m", "--meltBinRegion", dest="meltBinRegion",
        default="binC0304YCnt.merge/binC0304YCnt.merge.norm.chr_distrib.summary.xls", 
        help="File name for summary.xls output by meltRegionBinMatrix.py. One special symbol <C0304Y> should be put at the bin-size position. For example, file <bin10000Cnt.merge/bin10000Cnt.merge.norm.chr_distrib.summary.xls> should be <binC0304YCnt.merge/binC0304YCnt.merge.norm.chr_distrib.summary.xls> (which is default). The program will replace <C0304Y> with given bin sizes. Specifally giving <no> to ignore this part analysis.")
    parser.add_option("-b", "--bin-size", dest="bin_size",
        help="Bin size. Multiple sizes can be separated by < > or <,>.")
    parser.add_option("-t", "--tag", dest="tag",
        help="An unique string with only alohabets to label all labels in output Rmd.")
    #parser.add_option("-o", "--output-prefix", dest="out_prefix",
    #    help="The prefix of output files.")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='4_genome_bin_difference', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '4_genome_bin_difference'.")
    #parser.add_option("-d", "--doc-only", dest="doc_only",
    #    default=False, action="store_true", help="Specify to only generate doc.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.bin_size != None, "A string needed for -b"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    compare_summary = options.compare_summary
    meltBinRegion   = options.meltBinRegion
    bin_sizeL = re.split(r'[, ]*', options.bin_size)
    tag = options.tag
    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)

    print "# 样品比较 {#sample-compare-%s}\n" % tag
    print """
把基因组划分为相同大小的区域 (bin) (bin_size = ({})), 计算每个BIN的reads覆盖, 随后比较样品差异。结果如下。如果选取了不同大小的bin，每个bin分别呈现。    
""".format(', '.join(bin_sizeL))
    for bin in bin_sizeL:
        print '\n## Bin genome (bin_size={})\n'.format(bin)
        
        
        if meltBinRegion and meltBinRegion != 'no':

            print '\n### Chromosome reads distribution comparison for all samples (bin_size={})\n'.format(bin)

            curation_label = "Bin_chromsome_size_{}_read_distrib".format(bin)
            knitr_read_txt(report_dir, curation_label)
            tmp_meltBinRegion = meltBinRegion.replace('C0304Y', bin)
            pic_label_base = 'Bin-size-{}-chr-distrib-{}'.format(bin, tag)
            print """
绘制不同染色体上reads的分布，研究整体水平上样品间的reads丰度的异同。横轴代表染色体区域 (5'-3'), 纵轴代表每个bin的覆盖丰度 (CPM: reads count per million).

"""
            count = 0
            for line in open(tmp_meltBinRegion):
                count += 1
                png, chr = line.split()
                copypng(dest_dir, png)
                png = getRelativeDir(png, report_sub_dir)
                #png = pdf.replace('pdf', 'png')
                pic_label = pic_label_base + str(count)
                print '''

(ref:{pic_label}) Reads在{chr}染色体的分布 (bin_size={bin})。

```{{r {pic_label}, fig.cap="(ref:{pic_label})"}}
knitr::include_graphics("{png}")
```
'''.format(pic_label=pic_label, chr=chr, bin=bin, png=png)
        
        if compare_summary:
            print '\n### Differentiall detected regions (bin_size={})\n'.format(bin)
            curation_label = "Bin_chromsome_size_{}_sample_compare".format(bin)
            knitr_read_txt(report_dir, curation_label)
            
            pic_label_base = 'Bin-size-{}-samp-cmp-{}'.format(bin, tag)
            tmp_compare_summary = compare_summary.replace('C0304Y', bin)

            sumL = json_load(open(tmp_compare_summary))
            
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
根据基因组每个区间检测到的reads的丰度，从整体上对样品的相关性进行评估和样品分类，
确定相似性较高或较低的样品群 (Figure \@ref(fig:{pic_label_base}-pearson) and \@ref(fig:{pic_label_base}-pca))。

(ref:{pic_label_base}-pearson) 样品间Pearson相关性分析。图中每个格子代表纵轴和横轴两个样品的Pearson correlation coefficient (PCC). 格子的颜色表示了PCC值的大小。[pdf]({pdf}) [XLS]({xls})

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
样品间差异区域的鉴定和丰度展示 (FDR (false discovery rate) = {fdr}, log2FC (log2 fold change) = {log2fc})。

样品间差异区域数目统计，每个方格代表横轴样品相比于纵轴样品丰度显著高的区域的数目 (Figure \@ref(fig:{pic_label_base}-de-count))。

(ref:{pic_label_base}-de-count) 样品间差异区域数目统计，每个方格代表横轴样品相比于纵轴样品丰度显著高的区域的数目。[PDF]({pdf}) [XLS]({xls})

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

样品间差异区域丰度图谱。每一行为一个区域，每一列为一个样品，每个小方格代表对应区域在对应样品的丰度，其值的大小根据颜色展示 (Figure \@ref(fig:{pic_label_base}-de-profile))。

(ref:{pic_label_base}-de-profile) 样品间差异区域丰度图谱。每一行为一个区域，每一列为一个样品，每个小方格代表对应区域在对应样品的丰度，其值的大小根据颜色展示。[PDF]({pdf}) [XLS]({xls})

```{{r {pic_label_base}-de-profile, fig.cap="(ref:{pic_label_base}-de-profile)"}}
knitr::include_graphics("{png}")
```
""".format(pic_label_base=pic_label_base, png=png, pdf=pdf, 
    xls=xls, fdr=fdr, log2fc=log2fc)

            # DE compare -- each sample pair

            # The second element of sumL is a dict containing sample-compare 
            # information in sub-dict.
            compD = sumL[1]
            
            for comp_pair, comp_subD in compD.items():
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

{condA}和{condB}检测到的区域Reads丰度统计和区域注释 ([点击下载]({all}))。如果对应区域落在基因启动子区域 (转录起始位点上游1 kb和下游500 bp区域), 则以基因的注释作为此区域注释。一个区域可能对应多个基因。

{condA}相比于{condB}显著高的区域有{higherThan_Cnt}个，注释见[表格]({higherThan}).

{condB}相比于{condA}显著高的区域有{lowerThan_Cnt}个，注释见[表格]({lowerThan}).

""".format(condA=condA, condB=condB, bin=bin, all=all, higherThan_Cnt=higherThan_Cnt, lowerThan_Cnt=lowerThan_Cnt, higherThan=higherThan, lowerThan=lowerThan)

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


