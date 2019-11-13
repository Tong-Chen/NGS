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
    This is designed to summarize resulrts output by `geneBody_coverage2.py`.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import re
from tools import *
import pandas as pd
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
        metavar="FILEIN", help="`,` or ` ` separated a list of files. *.Log.final.out generated by `STAR` during mapping")
    parser.add_option("-l", "--labels", dest="label",
        metavar="LABEL", help="`,` or ` ` separated a list of labels to label each file. It must have same order as files.")
    parser.add_option("-o", "--output-prefix", dest="out_prefix",
        help="The prefix of output files.")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='2_mapping_quality', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '2_mapping_quality'.")
    parser.add_option("-d", "--doc-only", dest="doc_only",
        default=False, action="store_true", help="Specify to only generate doc.")
    parser.add_option("-n", "--number", dest="number", type="int", 
        default=40, help="Set the maximum allowed samples for barplot. Default 40.\
 If more than this number of samples are given, heatmap will be used instead. UNUSED.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readTwoColumnFile(fileL, labelL, header=0, index_col=0):
    tmpL = []
    for file, label in zip(fileL, labelL):
        print >>sys.stderr, file
        coverageM = pd.read_table(file, header=header, index_col=index_col)
        coverageM.columns = [label]
        tmpL.append(coverageM)
    mergeM = pd.concat(tmpL, axis=1)
    return mergeM
#-, ------------------------
def plot(fileL):
    for file in fileL:
        cmd = "s-plot barPlot -f " + file
        os.system(cmd)
#--------------------------------------
def plot_melt(total_melt, nameL):
    x_level = ["'"+i+"'" for i in nameL]
    x_level = '"'+','.join(x_level)+'"'
    cmd = ["s-plot barPlot -m TRUE -a Sample -R 90 -B set -O 1 -w 20 -u 25 -f ", 
            total_melt, ' -k free_y -L', x_level, 
            ' -y  \'Reads count or relative percent\' -x \'Samples\' ']
    #print ' '.join(cmd)
    os.system(' '.join(cmd))
#--------------------------------------
def plot_heatmap(totalTable, sample_count):
    height = sample_count // 2 + 1
    if height < 10:
        height = 10
    cmd = ["s-plot heatmapS -a TRUE -b TRUE -R TRUE", 
            "-x white -y blue -u 17 -F 12 -v", str(height), 
            "-f ", totalTable, "-I \"RPM\" -l top -T 2 -o log2 ", 
            "-Q \"c(10, 50, 90)\"", 
            "-S \"c('Gene start (first 20%)','Gene_middle (middle 80%)','Gene_end (last 20%)')\""]
    os.system(' '.join(cmd))

#---------------------------------------

def generateDoc(report_dir, report_sub_dir, totalTable, curation_label):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)
    pdf = totalTable+'.heatmapS.log2.pdf'
    copypdf(dest_dir, pdf)
    copy(dest_dir, totalTable)

    print "\n## 基因区Reads分布 {#Reads-distribution-along-gene-body}\n"
    curation_label = "Reads_distribution_along_gene_body"
    knitr_read_txt(report_dir,  curation_label)

    pdf = report_sub_dir+'/'+os.path.split(pdf)[-1]
    newTable = report_sub_dir+'/'+os.path.split(totalTable)[-1]
    png = pdf.replace('pdf', 'png')
    
    print """每个样品测序reads在基因区的分布如图 (Figure \@ref(fig:read-gene-body-fig))。用于评估建库mRNA的降解程度和用于测序reads的均衡性。

计算方式是把基因从5'到3'分为100个区间，计算每个区间的reads数。理想的情况是Reads在整个基因区从5'到3'分布一致，但通常情况下，基因两端的reads分布会少于中间部位，只要偏差不太大就可以继续使用。如果偏差很大，则应检查样品。

(ref:read-gene-body-fig) Summary average reads distribution along gene bodies. 每个样品对应的行颜色越均一、越深说明样品的测序质量越合格。 [PDF]({}) [XLS]({})\n
""".format(pdf, newTable)
    
    print '''```{{r read-gene-body-fig, fig.cap="(ref:read-gene-body-fig)"}}
knitr::include_graphics("{png}")
```
'''.format(png=png)

#--------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    fileL = re.split(r'[, ]*', file.strip())
    sample_readin = len(fileL)
    label = options.label
    labelL = re.split(r'[, ]*', label.strip())
    verbose = options.verbose
    op = options.out_prefix
    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    global debug
    debug = options.debug
    doc_only = options.doc_only
    num_samples_each_grp = options.number
    melt = 0
    if sample_readin <= num_samples_each_grp:
        melt = 1
    #-----------------------------------
    aDict = {}
    totalTable = op+".geneBody_coverage.xls"
    curation_label = os.path.split(sys.argv[0])[-1].replace('.', '_')
    if doc_only:
        generateDoc(report_dir, report_sub_dir, totalTable, curation_label)
        return 0
    mergeM = readTwoColumnFile(fileL, labelL, header=0, index_col=0)
    newIndex = [int(i)+1 for i in mergeM.index]
    mergeM.index = newIndex
    mergeM = mergeM.T   
    mergeM.to_csv(totalTable, sep=b"\t")


    plot_heatmap(totalTable, len(fileL))

    generateDoc(report_dir, report_sub_dir, totalTable, curation_label)
    ###--------multi-process------------------

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


