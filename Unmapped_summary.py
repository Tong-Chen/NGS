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
    This is designed to summarize results output by `parse_why_unmapp_blast.py`.
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
        metavar="FILEIN", help="`,` or ` ` separated a list of files. *.unmapped.table.parse.xls.dodgeBars.pdf in folder `unmapped`.")
    parser.add_option("-l", "--labels", dest="label",
        metavar="LABEL", help="`,` or ` ` separated a list of labels to label each file. It must have same order as files.")
    parser.add_option("-o", "--output-prefix", dest="out_prefix",
        help="The prefix of output files. UNUSED")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='2_mapping_quality', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '2_mapping_quality'.")
    parser.add_option("-d", "--doc-only", dest="doc_only",
        default=False, action="store_true", help="Specify to only generate doc. UNUSED.")
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
def plot_heatmap(totalTable):
    cmd = ["s-plot heatmapS -a TRUE -b TRUE -R TRUE", 
            "-x white -y blue -u 18 -v 30 -F 12 ", 
            "-f ", totalTable, "-I RPM"]
    os.system(' '.join(cmd))

#---------------------------------------

def generateDoc(report_dir, report_sub_dir, fileL, labelL, cntD):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)

    print "\n## 未比对reads的来源探索 {#Source-unmap-reads}\n"

    curation_label = "Unmapped_reads_exploring"
    knitr_read_txt(report_dir,  curation_label)
    
    print """
选取1000条未比对回基因组的reads，与NCBI的nt库进行比对，查看reads是否比对到其它物种或是测序错误造成的错配太多。每个reads在NCBI的nt库中选择一个最佳的匹配，若此最佳匹配能覆盖reads的50%以上，并且一致性达到80%，则认为此reads可能来源于这个匹配，从而确定reads的来源；反之，则认为reads未匹配到任何物种，可能是测序过程中引入的无关序列或测序错误。

"""
    
    len_fileL = len(fileL)
    group = 3
    
    for i in range(0, len_fileL, 3):
        pdfL = fileL[i:i+3]
        subL = labelL[i:i+3]
        copypdf(dest_dir, *pdfL)
        len_subF = len(pdfL)
        pdfL = [report_sub_dir+'/'+os.path.split(j)[-1] for j in pdfL]
        pngL = [j.replace('pdf', 'png') for j in pdfL]

        pdf_link = []  #[label_pdf](pdf), [label_pdf](pdf)
        for pdf, label in zip(pdfL, subL):
            tmp_155 = '['+label+'_pdf]'+'('+pdf+')'
            pdf_link.append(tmp_155)

        pdf_link = ' '.join(pdf_link)

        print "(ref:unmapp-origin-fig-{}) 未比对Reads的来源分布。From left to right, the samples are **{}**。{}\n".format(i, ', '.join(subL), pdf_link)
        
        pngFileL = []  #"png1", "png2", "png3"
        for png in pngL:
            tmp_164 = "'"+png+"'"
            pngFileL.append(tmp_164)
        pngFileL = ', '.join(pngFileL)
        print '''```{{r unmapp-origin-fig-{label}, out.width="{width}%", fig.cap="(ref:unmapp-origin-fig-{label})"}}
knitr::include_graphics(c({png}))
```
'''.format(label=i, png=pngFileL, width=int(100/len_subF))

#--------------------------------

def read_cnt_file(cnt_file):
    cntD = {}
    header = 1
    for line in open(cnt_file):
        if header:
            header -= 1
            continue
        lineL = line.strip().split('\t')
        sample = lineL[0]
        reads_cnt = lineL[-1]
        assert sample not in cntD, "Duplicate "+sample
        cntD[sample] = reads_cnt
    return cntD
#----------------------------

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
    #cnt_file = options.cnt
    # No use, can be deleted
    cntD = {}
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
    #curation_label = os.path.split(sys.argv[0])[-1].replace('.', '_')
    if doc_only:
        generateDoc(report_dir, report_sub_dir, fileL, labelL, cntD)
        return 0

    generateDoc(report_dir, report_sub_dir, fileL, labelL, cntD)
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


