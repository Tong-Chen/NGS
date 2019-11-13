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
    
    The program will try to identify related files based on <output_dir>, <typeL>, <label>.
    Generated files include <output_dir>/<type>.<label>.profile.pdf, 
                            <output_dir>/<type>.<profile.gene>.profile.kmeans.pdf,
                            <output_dir>/<type>.<profile.gene>.heatmap.kmeans.pdf.

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
    parser.add_option("-p", "--input-dir", dest="output_dir",
        help="Normally $(summary_dir)/")
    parser.add_option("-t", "--typeL", dest="typeL",
        help="`,` or ` ` separated strings.")
    parser.add_option("-l", "--label", dest="label",
        default='profile.gene', 
        help="A string to represent the result. Default <profile.gene>.")
    #parser.add_option("-o", "--output-file-output_dir", dest="op",
    #    help="Specify output file output_dir")
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
    assert options.output_dir != None, "A filename needed for -p"
    return (options, args)
#--------------------------------------------------------------------


def generateDoc(report_dir, report_sub_dir, output_dir, tag, typeL):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)
    
    print "\n## 基因区结合图谱绘制  {#ChIP-sample-comp-gene-profile}\n"
    
    curation_label = "Gene_region_line_heatmap"
    knitr_read_txt(report_dir,  curation_label)
    
    print """

取基因转录起始区域上游3 Kb，转录终止区域下游 3 Kb，分为大小为100 bp的bin。同时把基因区根据长度统一分为50个bin，计算每个bin区域富集的reads相对于input的reads的log2 ratio，以线图和热图的方式展示，从而方便比较样品内和样品间的结合图谱。


"""
    count = 0
    for type in typeL:
        count += 1
        label = tag.replace('.', '') + str(count)


        profile_pdf = output_dir + '/' + type + '.' +tag + '.profile.pdf'
        profile_kmeans_pdf = output_dir + '/' + type + '.' +tag + '.profile.kmeans.pdf'
        heatmap_kmeans_pdf = output_dir + '/' + type + '.' +tag + '.heatmap.kmeans.pdf'

        copypdf(dest_dir, profile_pdf, profile_kmeans_pdf, heatmap_kmeans_pdf)

        profile_pdf,profile_kmeans_pdf,heatmap_kmeans_pdf = getRelativeDir([profile_pdf,profile_kmeans_pdf,heatmap_kmeans_pdf],report_sub_dir)
        profile_png = profile_pdf.replace('pdf', 'png')
        profile_kmeans_png = profile_kmeans_pdf.replace('pdf', 'png')
        heatmap_kmeans_png = heatmap_kmeans_pdf.replace('pdf', 'png')
        profile_pdf_n,profile_kmeans_pdf_n,heatmap_kmeans_pdf_n = getFileName([profile_pdf,profile_kmeans_pdf,heatmap_kmeans_pdf])


        print """
### {type} 结合图谱

基因区总的结合图谱如 Figure \@ref(fig:{label}-profile-{count})。

(ref:{label}-profile-{count}) Binding profile along upstream and downstream 3 kb regions of gene body as well as gene body regions. [{profile_pdf_n}]({profile_pdf}) 

```{{r {label}-profile-{count}, fig.cap="(ref:{label}-profile-{count})"}}
knitr::include_graphics("{profile_png}")
```

根据基因区结合图谱分布模式，聚类为2类的图谱如 Figure \@ref(fig:{label}-profile-{count}-kmeans)。

(ref:{label}-profile-{count}-kmeans) Two clusters of binding profiles along upstream and downstream 3 kb regions of gene body as well as gene body regions. The clusters were generated using k-means on binding profiles. [{profile_kmeans_pdf_n}]({profile_kmeans_pdf}) 

```{{r {label}-profile-{count}-kmeans, fig.cap="(ref:{label}-profile-{count}-kmeans)"}}
knitr::include_graphics("{profile_kmeans_png}")
```

根据基因区结合图谱分布模式，聚类为4类的Heatmap图谱如 Figure \@ref(fig:{label}-heatmap-{count}-kmeans)。

(ref:{label}-heatmap-{count}-kmeans) Four clusters of binding profiles along upstream and downstream 3 kb regions of gene body as well as gene body regions. The clusters were generated using k-means on binding patterns. [{heatmap_kmeans_pdf_n}]({heatmap_kmeans_pdf}) 

```{{r {label}-heatmap-{count}-kmeans, fig.cap="(ref:{label}-heatmap-{count}-kmeans)"}}
knitr::include_graphics("{heatmap_kmeans_png}")
```

""".format(label=label, type=type, count=count,  
       profile_pdf_n=profile_pdf_n, profile_pdf=profile_pdf, profile_png=profile_png, 
       profile_kmeans_pdf_n=profile_kmeans_pdf_n, profile_kmeans_pdf=profile_kmeans_pdf, profile_kmeans_png=profile_kmeans_png, 
       heatmap_kmeans_pdf_n=heatmap_kmeans_pdf_n, heatmap_kmeans_pdf=heatmap_kmeans_pdf, heatmap_kmeans_png=heatmap_kmeans_png) 


#------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    output_dir = options.output_dir
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

    generateDoc(report_dir, report_sub_dir, output_dir, label, typeL) 

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


