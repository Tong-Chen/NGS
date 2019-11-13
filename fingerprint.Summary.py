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
    This is designed to summary the file output of plotFingerprint.
    
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
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="<,> or < > separated a list of files.")
    parser.add_option("-l", "--labels", dest="label",
        metavar="LABEL", help="`,` or ` ` separated a list of labels to label each file. It must have same order as files.")
    #parser.add_option("-b", "--bed", dest="bed",
    #    help="The bed file given to <bedtools coverage -a > ('main_chrom_bed' in this example). The forth coulumn will be used as region names.")
    #parser.add_option("-o", "--output-file-prefix", dest="op",
    #    help="Specify output file prefix")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='3_ChIP_quality', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '3_ChIP_quality'.")
    #parser.add_option("-a", "--appendFile", dest="append",
    #    default='', help="A list of files to be appended. Optional.")
    parser.add_option("-d", "--doc-only", dest="doc_only",
        default=False, action="store_true", help="Specify to only generate doc.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def generateDoc(report_dir, report_sub_dir, pdfL, labelL):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)
    images = "images"
    os.system('mkdir -p '+report_dir+'/'+images)
    copypng(report_dir+'/'+images, "/MPATHB/self/resource/sample/QC_fingerprint.png")
    
    #copy(dest_dir, *fileL)
    copypdf(dest_dir, *pdfL)
    
    pdfL = getRelativeDir(pdfL, report_sub_dir)
    pngL = [pdf.replace('pdf', 'png') for pdf in pdfL]

    print "\n## 抗体富集效率评估 {#ChIP-quality-IP}\n"
    
    curation_label = "ChIP_enrichment_efficiency"
    knitr_read_txt(report_dir,  curation_label)

    print """
采用FingerPrint展示ChIP富集的信号与背景信号的区分度 (Figure \@ref(fig:chip-ip-sample))。绘制方法是，把基因组分成很小的区间，统计每个区间上reads的数目，并且按测序覆盖度由低到高排序，作为横轴；纵轴为这些区域的累计reads数占测序总reads数的比例。比如点 (0.6, 0.2)表示60%的基因组区域覆盖的reads数占总reads数的20%。
理想条件下，input样品的reads在基因组上是均匀分布的，绘制出来的应该是一条直线；如果对基因组的覆盖达到100%，应该是从原点处出发、斜率为1的对角直线。
一个非常特异并且富集程度特别高的ChIP实验中，测序得到的reads集中在基因组的少数区域，大部分区域没有覆盖，绘制出来的应该是平缓向右延伸，而后接近直角的突然上升曲线。

对于转录因子等结合强度高并且结合区域小的ChIP实验，富集质量高的样品FingerPrint图的区分度会很明显。而对组蛋白修饰等结合区域比较广的ChIP实验，则FingerPrint图的区分度不太明显。反过来，从FingerPrint图的分布，也可以判断ChIP信号的类型，sharp型或者broad型。

```{r chip-ip-sample, fig.cap="Examples to show you how the nature of the ChIP signal (narrow and high vs. wide and not extremely high) is reflected in the 'fingerprint' plots."}
knitr::include_graphics("images/QC_fingerprint.png")
```
"""
    count = 1
    for pdf, png, sample in zip(pdfL, pngL, labelL):
        label = "chip-ip-quality-"+str(count)
        count += 1
        print """
(ref:{label}) Fingerprints for sample {sample} to show IP-efficiency or signal type. [PDF]({pdf})

```{{r {label}, fig.cap="(ref:{label})"}}
knitr::include_graphics("{png}")
```

""".format(png=png, label=label, pdf=pdf, sample=sample)

#------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    fileL = re.split(r'[, ]*', file.strip())
    label = options.label
    labelL = re.split(r'[, ]*', label.strip())
    verbose = options.verbose
    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    doc_only = options.doc_only

    global debug
    debug = options.debug
    #-----------------------------------

    generateDoc(report_dir, report_sub_dir, fileL, labelL) 

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


