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
    This is designed to summary the file output by bedtools_genomecov.py.
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
    parser.add_option("-o", "--output-file-prefix", dest="op",
        help="Specify output file prefix")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='2_mapping_quality', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '2_mapping_quality'.")
    parser.add_option("-a", "--appendFile", dest="append",
        default='', help="A list of files to be appended. Optional.")
    parser.add_option("-d", "--doc-only", dest="doc_only",
        default=False, action="store_true", help="Specify to only generate doc.")
    parser.add_option("-t", "--test-line", dest="test",
        type="int", default=0, help="Give a positive number (n) to test the program using on the first n reads. Default 0 meaning using all reads in file.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readFile(file, label, op_fh):
    header = 1
    for line in open(file):
        if header:
            header -= 1
            continue
        lineL = line.strip().split('\t')
        for i in range(int(lineL[2])):
            print >>op_fh, '{}\t{}\t{}'.format(label, lineL[1], 'DNA fragment length')
        print >>op_fh, '{}\t{}\t{}'.format(label, lineL[2], 'DNA fragment count')
#---------------------------
def cntFile(file, test=0):
    header = 1
    lengthD = {}
    cntD = {}
    baseD = {}
    len_line = 0
    read_in_cnt = 0
    for line in open(file):
        lineL = line.strip().split('\t')
        read_in_cnt += 1
        if header:
            headerL = lineL
            len_line = len(lineL)
            for i in range(7, len_line):
                baseD[lineL[i]] = []
            header -= 1
            continue
        length   = int(lineL[3])
        cnt   = int(lineL[4])
        lengthD[length] = lengthD.get(length, 0)+cnt
        cntD[cnt] = cntD.get(cnt, 0)+1
        for i in range(7, len_line):
            for j in range(cnt):
                baseD[headerL[i]].append(lineL[i])
        if test and (read_in_cnt>test):
            break
    #--------------------------------------
    sum_reads = sum(lengthD.values())

    #assert sum_reads == sum(cntD.keys()), file
    lengthFreqD = {}
    for length,length_cnt in lengthD.items():
        lengthFreqD[length] = length_cnt / sum_reads * 100
    cntFreqD = {}
    sum_cnt = sum(cntD.values())
    test_sum = 0
    for cnt, cnt_cnt in cntD.items():
        test_sum += cnt * cnt_cnt
        cntFreqD[cnt] = cnt_cnt / sum_cnt * 100
    assert sum_reads == test_sum, file
    return lengthD, lengthFreqD, cntD, cntFreqD, baseD
#---------------------------------------------

#----------------------------------

def generateDoc(report_dir, report_sub_dir, cov_file, 
        curation_label, appendL, labelL):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)
    cov_file_pdf = cov_file + ".dodgeBars.pdf"
    if appendL:
        copy(dest_dir, *appendL)

    copy(dest_dir, cov_file)
    copypdf(dest_dir, cov_file_pdf)

    print "\n## 基因组测序覆盖度和测序深度评估 {#seq-breadth-depth}\n"
    
    curation_label = "genome_coverage_depth"
    knitr_read_txt(report_dir,  curation_label)

    print """
**Reads sheathed genome regions (%)**: 测到的基因组区域占总基因组区域的比例，测序覆盖度。

**Sequencing depth relative to whole genome**: 测序reads相对于全基因组区域的覆盖度，测序深度。

**Sequencing depth relative to sheathed genome**: 测序reads相对于测到的基因组区域的覆盖度，测序深度。

"""

    cov_file_pdf = report_sub_dir+'/'+os.path.split(cov_file_pdf)[-1]
    cov_file = report_sub_dir+'/'+os.path.split(cov_file)[-1]
    cov_file_png = cov_file_pdf.replace('pdf', 'png')

    print "(ref:cov-file-distrib-fig) 测序覆盖度和测序深度评估。**Reads sheathed genome regions (%)**: 测到的基因组区域占总基因组区域的比例；**Sequencing depth relative to whole genome**: 测序reads相对于全基因组区域的覆盖度，测序深度；**Sequencing depth relative to sheathed genome**: 测序reads相对于测到的基因组区域的覆盖度，测序深度。[PDF]({}) [XLS]({})\n".format(cov_file_pdf, cov_file)

    print """```{{r cov-file-distrib-fig, fig.cap="(ref:cov-file-distrib-fig)"}}
knitr::include_graphics("{png}")
```
""".format(png=cov_file_png)


    #if appendL:
    #    appendList = generateLink(appendL, labelL, 'seq_fragment_attribute', 
    #            report_sub_dir, '\n* ')
    #    print """## 测序片段统计表格下载 {{#seq-frag-xls}}

#* {}
#""".format(appendList)
#-------------------------------

def output(fileL, fh):
    for file in fileL:
        header = 1
        for line in open(file):
            if header:
                header -= 1
                continue
            print >>fh, line,
#---------------------------------------------------------------
def plot(cov_file, nameL):
    x_level = ["'"+i+"'" for i in nameL]
    x_level = '"'+','.join(x_level)+'"'
    cmd = ['s-plot barPlot -m TRUE -a sample -d dodge -P none', 
           '-B variable -O 1 -k free_y -R 45', '-f', cov_file, 
           '-L', x_level]
    print >>sys.stderr, ' '.join(cmd)
    os.system(' '.join(cmd))
#---------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    fileL = re.split(r'[, ]*', file.strip())
    label = options.label
    labelL = re.split(r'[, ]*', label.strip())
    verbose = options.verbose
    test  = options.test
    op = options.op
    append  = options.append
    if append:
        appendL = re.split(r'[, ]*', append.strip())
    else:
        appendL = []
    cov_file = op+'.genomecov_summary.xls'

    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    doc_only = options.doc_only
    curation_label = os.path.split(sys.argv[0])[-1].replace('.',  '_')
    if doc_only:
        generateDoc(report_dir, report_sub_dir, cov_file, 
                curation_label, appendL, labelL)
        return 0
    global debug
    debug = options.debug
    #-----------------------------------

    cov_file_fh = open(cov_file, 'w')

    print >>cov_file_fh, "variable\tvalue\tsample"
    
    output(fileL, cov_file_fh)

    cov_file_fh.close()

    plot(cov_file, labelL)

    generateDoc(report_dir, report_sub_dir, cov_file, 
            curation_label, appendL, labelL)

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


