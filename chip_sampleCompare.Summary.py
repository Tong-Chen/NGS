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
    This is designed to summary results of sample compare parson correlations and PCA results.
    
    The program will try to identify related files based on <prefix>, <label>, <typeL>.
    Generated files include <prefix>.<labeL>.correlation.xls, 
                            <prefix>.<labeL>.correlation.pdf,
                            <prefix>.<label>.xls.pca.scale.pdf, 
                            <prefix>.<type1>.<labeL>.correlation.xls, 
                            <prefix>.<type1>.<labeL>.correlation.pdf,
                            <prefix>.<type1>.<label>.xls.pca.scale.pdf, 

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
    parser.add_option("-p", "--input-prefix", dest="prefix",
        help="Normally $(summary_dir)/$(prefix)")
    parser.add_option("-l", "--labels", dest="label",
        default="multiBigwigSummary", help="Compare type. Currently <multiBigwigSummary> (default) is supported. ")
    parser.add_option("-t", "--typeL", dest="typeL",
        help="`,` or ` ` separated strings.")
    #parser.add_option("-o", "--output-file-prefix", dest="op",
    #    help="Specify output file prefix")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='4_Sample_compare', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '4_Sample_compare'.")
    #parser.add_option("-a", "--appendFile", dest="append",
    #    default='', help="A list of files to be appended. Optional.")
    parser.add_option("-d", "--doc-only", dest="doc_only",
        default=False, action="store_true", help="Specify to only generate doc.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.prefix != None, "A filename needed for -p"
    return (options, args)
#--------------------------------------------------------------------


def generateDoc(report_dir, report_sub_dir, prefix, tag, typeL):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)
    
    
    #copy(dest_dir, *fileL)
    #copypdf(dest_dir, *pdfL)
    
    print "\n## 样品整体比较和分类  {#ChIP-sample-comp-pca-cor}\n"
    
    curation_label = "Sample_compare_classification"
    knitr_read_txt(report_dir,  curation_label)
    
    print """

把基因组分为大小为1 kb的bin，统计每个bin的标准化后（如果有input，则为相对于input）的reads丰度，然后计算样品之间的Spearman相关系数和PCA主成分分聚类，判断样品之间的整体异同和分类趋势。
"""
    cor_xls = '.'.join([prefix, tag, "correlation.xls"]) 
    cor_pdf = '.'.join([prefix, tag, "correlation.pdf"])
    pca_pdf = '.'.join([prefix, tag, "xls.pca.scale.pdf"])

    copy(dest_dir, cor_xls)
    copypdf(dest_dir, cor_pdf, pca_pdf)

    cor_xls, cor_pdf, pca_pdf = getRelativeDir([cor_xls, cor_pdf,  pca_pdf], report_sub_dir)
    cor_png = cor_pdf.replace('pdf', 'png')
    pca_png = pca_pdf.replace('pdf', 'png')

    count = 1
    label = "ChIP-sample-comp-pca-cor-"+str(count)
    print """
(ref:{label}) Spearman correlation (left, [PDF]({cor_pdf}), [XLS]({cor_xls})) and sample classification (right, [PDF]({pca_pdf})) showing sample similarity of ChIP profiles genome-widely.

```{{r {label}, fig.cap="(ref:{label})"}}
knitr::include_graphics(c("{cor_png}", "{pca_png}"))
```

""".format(cor_png=cor_png, label=label, cor_pdf=cor_pdf, pca_png=pca_png, pca_pdf=pca_pdf, cor_xls=cor_xls)

    if not typeL:
        return 
    for type in typeL:
        count += 1
        label = "ChIP-sample-comp-pca-cor-"+str(count)

        cor_xls = '.'.join([prefix, type, tag, "correlation.xls"]) 
        cor_pdf = '.'.join([prefix, type, tag, "correlation.pdf"])
        pca_pdf = '.'.join([prefix, type, tag, "xls.pca.scale.pdf"])

        copy(dest_dir, cor_xls)
        copypdf(dest_dir, cor_pdf, pca_pdf)

        cor_xls, cor_pdf, pca_pdf = getRelativeDir([cor_xls, cor_pdf,  pca_pdf], report_sub_dir)
        cor_png = cor_pdf.replace('pdf', 'png')
        pca_png = pca_pdf.replace('pdf', 'png')

        print """
## 样品整体比较和分类 ({type})

(ref:{label}) Spearman correlation (left, [PDF]({cor_pdf}), [XLS]({cor_xls})) and sample classification (right, [PDF]({pca_pdf})) showing sample similarity of {type} ChIP profiles genome-widely.

```{{r {label}, fig.cap="(ref:{label})"}}
knitr::include_graphics(c("{cor_png}", "{pca_png}"))
```

""".format(cor_png=cor_png, label=label, cor_pdf=cor_pdf, pca_png=pca_png, pca_pdf=pca_pdf, type=type, cor_xls=cor_xls)


#------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    prefix = options.prefix
    label = options.label
    typeL = options.typeL.strip()
    if typeL:
        typeL = re.split(r'[, ]*', typeL)
    else:
        typeL = []
    verbose = options.verbose
    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    doc_only = options.doc_only

    global debug
    debug = options.debug
    #-----------------------------------

    generateDoc(report_dir, report_sub_dir, prefix, label, typeL) 

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


