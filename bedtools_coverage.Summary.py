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
    This is designed to summary the file output of 
    bedtools coverage -a $(main_chrom_bed) -b final.bam -bed
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
    parser.add_option("-b", "--bed", dest="bed",
        help="The bed file given to <bedtools coverage -a > ('main_chrom_bed' in this example). The forth coulumn will be used as region names.")
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
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def output(fileL, labelL, nameL, cov_file_raw_fh, cov_file_raw_normalize_fh, 
           cov_file_percent_fh, cov_file_percent_normalize_fh, 
           cov_file_coverage_fh):
    for file, label in zip(fileL, labelL):
        cntD = {}
        lengthD = {}
        covD = {}
        for line in open(file):
            lineL = line.split()
            key = lineL[3]
            assert key not in cntD, "Duplicate "+key+' for file '+file
            cntD[key] = int(lineL[4])
            lengthD[key] = int(lineL[6])
            covD[key] = lineL[7]
        cntL = [cntD.get(name, 0) for name in nameL]
        print >>cov_file_raw_fh, \
            "{}\t{}".format(label, '\t'.join([str(i) for i in cntL]))
        cnt_normalize = [cntD.get(name, 0)/lengthD.get(name, 1) for name in nameL]
        print >>cov_file_raw_normalize_fh, \
            "{}\t{}".format(label, '\t'.join([str(i) for i in cnt_normalize]))

        readsSum = sum(cntL)
        percentL = [100.0*i/readsSum for i in cntL]
        print >>cov_file_percent_fh, \
            "{}\t{}".format(label, '\t'.join([str(i) for i in percentL]))
        cnt_normalize_sum = sum(cnt_normalize)
        norm_percentL = [100*i/cnt_normalize_sum for i in cnt_normalize]
        print >>cov_file_percent_normalize_fh, \
            "{}\t{}".format(label, '\t'.join([str(i) for i in norm_percentL]))
        
        covL = [covD.get(name, 0) for name in nameL]
        print >>cov_file_coverage_fh, \
            "{}\t{}".format(label, '\t'.join([str(i) for i in covL]))

#-----------------------------------------------------------

def generateDoc(report_dir, report_sub_dir, curation_label, appendL, labelL, 
        cov_file_raw, cov_file_raw_normalize, cov_file_percent, 
        cov_file_percent_normalize, cov_file_coverage):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)
    
    pdfL = [cov_file_raw+'.heatmapS.log2.pdf', 
            cov_file_raw_normalize+'.heatmapS.log2.pdf', 
            cov_file_percent+'.heatmapS.pdf', 
            cov_file_percent_normalize+'.heatmapS.pdf', 
            cov_file_coverage+'.heatmapS.pdf'
            ]
    fileL = [cov_file_raw, cov_file_raw_normalize, cov_file_percent, 
        cov_file_percent_normalize, cov_file_coverage]

    labelL = ["cov-file-raw-fig", "cov-file-raw-normalize-fig", 
              "cov-file-percent-fig", 
              "cov-file-percent-normalize-fig", "cov-file-coverage-fig"]
    count = len(fileL)

    if appendL:
        copy(dest_dir, *appendL)

    copy(dest_dir, *fileL)
    copypdf(dest_dir, *pdfL)

    print "\n## 染色体测序覆盖度和测序深度评估 {#chromosome-seq-breadth-depth}\n"
    
    curation_label = "chromosome_coverage_depth"
    knitr_read_txt(report_dir,  curation_label)

    print """
为了评估不同染色体的reads分布，统计了每个染色体比对到的reads的绝对数目 (Figure \@ref(fig:{}))和reads的相对比例 (Figure \@ref(fig:{}))以及测到的染色体占总染色体的比例 (染色体的覆盖度) (Figure \@ref(fig:{}))。同时考虑到不同染色体长度不一，长度越长的染色体被测到的机会也越大，因此根据染色体的长度对比对到每个染色体的reads的绝对数目和reads的相对比例做了标准化 (normalize) (Figure \@ref(fig:{}) and \@ref(fig:{}))。

""".format(labelL[0], labelL[2], labelL[4], labelL[1], labelL[3])

    annoL= ["比对到每条染色体的reads的绝对数目统计。", 
            "比对到每条染色体的reads的绝对数目统计 (染色体长度标准化)。", 
            "比对到每条染色体的reads的相对比例。", 
            "比对到每条染色体的reads的相对比例 (染色体长度标准化)。", 
            "染色体覆盖度。"]
    
    for i in range(count):
        cov_file_pdf = getRelativeDir(pdfL[i], report_sub_dir)
        cov_file = getRelativeDir(fileL[i], report_sub_dir)
        cov_file_png = cov_file_pdf.replace('pdf', 'png')

        print "(ref:{}) {} [PDF]({}) [XLS]({})\n"\
                .format(labelL[i], annoL[i], cov_file_pdf, cov_file)

        print """```{{r {label}, fig.cap="(ref:{label})"}}
knitr::include_graphics("{png}")
```
""".format(png=cov_file_png, label=labelL[i])


    #if appendL:
    #    appendList = generateLink(appendL, labelL, 'seq_fragment_attribute', 
    #            report_sub_dir, '\n* ')
    #    print """## 测序片段统计表格下载 {{#seq-frag-xls}}

#* {}
#""".format(appendList)
#-------------------------------

def plot(sample_cnt, cov_file_raw, cov_file_raw_normalize, 
           cov_file_percent, cov_file_percent_normalize, 
           cov_file_coverage):
    height = sample_cnt // 2 + 2
    if height < 10:
        height = 10
    cmd = ["s-plot heatmapS -a TRUE -A 45 -b TRUE -R TRUE", 
            "-x white -y blue -u 22 -F11 -T 2 -l top -o log2", 
            "-f ", cov_file_raw, "-I \"Count\"", '-v', str(height), 
            '-t "Number of reads mapped to each chromosome"']
    os.system(' '.join(cmd))
    cmd = ["s-plot heatmapS -a TRUE -A 45 -b TRUE -R TRUE", 
            "-x white -y blue -u 22 -F11 -T 2 -l top -o log2", 
            "-f ", cov_file_raw_normalize, "-I \"Count (norm)\"", 
            '-v', str(height), 
            '-t "Number of reads mapped to each chromosome (normalized by chromosome length)"']
    os.system(' '.join(cmd))
    cmd = ["s-plot heatmapS -a TRUE -A 45 -b TRUE -R TRUE", 
            "-x white -y blue -u 22 -F11 -T 2 -l top", 
            "-f ", cov_file_percent, "-I \"Percent (%)\"", '-v', str(height), 
            '-t "Percent of reads mapped to each chromosome"']
    os.system(' '.join(cmd))
    cmd = ["s-plot heatmapS -a TRUE -A 45 -b TRUE -R TRUE", 
            "-x white -y blue -u 22 -F11 -T 2 -l top", 
            "-f ", cov_file_percent_normalize, "-I \"Percent (norm) (%)\"", '-v', str(height), 
            '-t "Percent of reads mapped to each chromosome (normalized by chromosome length)"']
    os.system(' '.join(cmd))
    cmd = ["s-plot heatmapS -a TRUE -A 45 -b TRUE -R TRUE", 
            "-x white -y blue -u 22 -F11 -T 2 -l top", 
            "-f ", cov_file_coverage, "-I 'Coverage (%)'", '-v', str(height), 
            '-t "Percent of sheathed chromosomes"']
    os.system(' '.join(cmd))
    
#-----------------------------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    fileL = re.split(r'[, ]*', file.strip())
    label = options.label
    labelL = re.split(r'[, ]*', label.strip())
    bed = options.bed
    verbose = options.verbose
    op = options.op
    append  = options.append
    if append:
        appendL = re.split(r'[, ]*', append.strip())
    else:
        appendL = []
    cov_file_raw = op+'.chr_cov_raw.xls'
    cov_file_raw_normalize = op+'.chr_cov_raw_normalize.xls'
    cov_file_percent = op+'.chr_cov_percent.xls'
    cov_file_percent_normalize = op+'.chr_cov_percent_normalize.xls'
    cov_file_coverage = op + '.chr_cov_region.xls'

    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    doc_only = options.doc_only
    curation_label = os.path.split(sys.argv[0])[-1].replace('.',  '_')
    if doc_only:
        generateDoc(report_dir, report_sub_dir, curation_label, appendL, labelL, 
            cov_file_raw, cov_file_raw_normalize, cov_file_percent, 
            cov_file_percent_normalize, cov_file_coverage)
        return 0
    global debug
    debug = options.debug
    #-----------------------------------

    cov_file_raw_fh = open(cov_file_raw, 'w')
    cov_file_raw_normalize_fh = open(cov_file_raw_normalize, 'w')
    cov_file_percent_fh = open(cov_file_percent, 'w')
    cov_file_percent_normalize_fh = open(cov_file_percent_normalize, 'w')
    cov_file_coverage_fh = open(cov_file_coverage, 'w')
    
    nameL = [line.split()[3] for line in open(bed)]

    print >>cov_file_raw_fh, "Sample\t{}".format('\t'.join(nameL))
    print >>cov_file_raw_normalize_fh, "Sample\t{}".format('\t'.join(nameL))
    print >>cov_file_percent_fh, "Sample\t{}".format('\t'.join(nameL))
    print >>cov_file_percent_normalize_fh, "Sample\t{}".format('\t'.join(nameL))
    print >>cov_file_coverage_fh, "Sample\t{}".format('\t'.join(nameL))
    
    output(fileL, labelL, nameL, cov_file_raw_fh, cov_file_raw_normalize_fh, 
           cov_file_percent_fh, cov_file_percent_normalize_fh, 
           cov_file_coverage_fh)

    cov_file_raw_fh.close()
    cov_file_raw_normalize_fh.close()
    cov_file_percent_fh.close()
    cov_file_percent_normalize_fh.close()
    cov_file_coverage_fh.close()

    plot(len(fileL), cov_file_raw, cov_file_raw_normalize, 
           cov_file_percent, cov_file_percent_normalize, 
           cov_file_coverage)

    generateDoc(report_dir, report_sub_dir, curation_label, appendL, labelL, 
            cov_file_raw, cov_file_raw_normalize, cov_file_percent, 
            cov_file_percent_normalize, cov_file_coverage)

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


