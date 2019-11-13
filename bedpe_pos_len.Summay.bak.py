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
    This is designed to summary the file output by bedpe_pos_len.py.
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
        default='3_read_length_distrib', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '2_mapping_quality'.")
    parser.add_option("-a", "--appendFile", dest="append",
        default='', help="A list of files to be appended.")
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

def plot(op):
    cmd = ['s-plot densityHistPlot -m TRUE -s set -S 1 -D free -v TRUE -f',
            op, '-y Frequency -x "DNA fragment length (bp) or count"', 
            '-u 30']
    os.system(' '.join(cmd))
#----------------------------------

def generateDoc(report_dir, report_sub_dir, frag_len, frag_cnt, frag_base, 
        curation_label, appendL, labelL):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)
    frag_len_pdf = frag_len + ".lines.pdf"
    frag_len_smooth_pdf = frag_len + ".lines.smooth.pdf"
    frag_cnt_pdf = frag_cnt + ".lines.pdf"
    frag_base_pdf = frag_base+'.boxplot.violin_nb.pdf'
    copy(dest_dir, *appendL)
    copypdf(dest_dir, frag_len_pdf, frag_len_smooth_pdf, frag_cnt_pdf, frag_base_pdf)

    print "\n## Length and count distribution of sequenced fragments\n"
    
    knitr_read_txt(report_dir,  curation_label)

    print """测序片段长度和数目分布。

为了保证结果的准确性，只选取DNA双端测序Reads都正确的比对到基因组唯一位置的片段做后续分析。片段的长度定义为双端Reads所覆盖的基因组区域。

Reads mapped to a unique genomic location were used for downstream analyses. Paired-end reads aligned to the same chromosome with a correct orientation were retained for downstream size analyses.

1. DNA测序片段大小分布如 Figure \@ref(fig:frag-len-distrib-fig)所示；

   * **Fragment length distrib (count)** 表示不同大小DNA片段的绝对数目的分布，不同样品间会受到测序丰度的影响；横轴表示片段的大小，纵轴表示特定大小片段测到的次数。举个例子，点(150, 200)表示长度为150的片段的出现次数为200次。
   * **Fragment length distrib (frequency (%))** 表示不同大小DNA片段的频率分布的分布；横轴表示片段的大小，纵轴表示特定大小片段测到的相对频率。举个例子，点(150, 2)表示长度为150的片段的出现频率为2%。
"""

    frag_len_pdf = report_sub_dir+'/'+os.path.split(frag_len_pdf)[-1]
    frag_len_png = frag_len_pdf.replace('pdf', 'png')

    print "(ref:frag-len-distrib-fig) DNA测序片段大小分布。**Fragment length distrib (count)** 表示不同大小DNA片段的绝对数目的分布，不同样品间会受到测序丰度的影响；**Fragment length distrib (frequency (%))** 表示不同大小DNA片段的频率分布的分布。[PDF]({})\n".format(frag_len_pdf)

    print """```{{r frag-len-distrib-fig, fig.cap="(ref:frag-len-distrib-fig)"}}
knitr::include_graphics("{png}")
```
""".format(png=frag_len_png)

    print """ 
2. 线性拟合^[线性拟合是根据数据的分布特征拟合出一条与其变化趋势相近的平滑曲线，会使展示结果更清晰，但也会略掉部分细节信息。]的DNA测序片段大小分布如 Figure \@ref(fig:frag-len-smooth-distrib-fig)所示；

   * **Fragment length distrib (count)** 表示不同大小DNA片段的绝对数目的分布，不同样品间会受到测序丰度的影响；
   * **Fragment length distrib (frequency (%))** 表示不同大小DNA片段的频率分布的分布； 
"""

    frag_len_smooth_pdf = report_sub_dir+'/'+os.path.split(frag_len_smooth_pdf)[-1]
    frag_len_smooth_png = frag_len_smooth_pdf.replace('pdf', 'png')

    print "(ref:frag-len-smooth-distrib-fig) DNA测序片段大小分布 (拟合结果)。**Fragment length distrib (count)** 表示不同大小DNA片段的绝对数目的分布，不同样品间会受到测序丰度的影响；**Fragment length distrib (frequency (%))** 表示不同大小DNA片段的频率分布的分布。[PDF]({})\n".format(frag_len_smooth_pdf)

    print """```{{r frag-len-smooth-distrib-fig, fig.cap="(ref:frag-len-smooth-distrib-fig)"}}
knitr::include_graphics("{png}")
```
""".format(png=frag_len_smooth_png)
   
    print """ 
3. DNA序列片段测到的频率分布如 Figure \@ref(fig:frag-cnt-distrib-fig)所示：

   * **Fragment occurence number distrib (count)** 表示DNA片段被测到的次数的分布。横轴为片段的出现次数，纵轴为特定次数的片段的条数。比如点 (100, 5)表示有5条片段各被测到100次。

   * **Fragment occurence number distrib (frequency (%))** 表示DNA片段被测到的次数的分不。横轴为片段的出现次数，纵轴为特定次数的片段的相对比例。比如点 (100, 20)表示有20%的片段被测到100次。

"""

    frag_cnt_pdf = report_sub_dir+'/'+os.path.split(frag_cnt_pdf)[-1]
    frag_cnt_png = frag_cnt_pdf.replace('pdf', 'png')

    print "(ref:frag-cnt-distrib-fig) 测序得到的DNA序列片段的频率分布。**Fragment occurence number distrib (count)** 表示DNA片段被测到的次数的分布。横轴为片段的出现次数，纵轴为特定次数的片段的条数; **Fragment occurence number distrib (frequency (%))** 表示DNA片段被测到的次数的分不。横轴为片段的出现次数，纵轴为特定次数的片段的相对比例。 [PDF]({})\n".format(frag_cnt_pdf)

    print """```{{r frag-cnt-distrib-fig, fig.cap="(ref:frag-cnt-distrib-fig)"}}
knitr::include_graphics("{png}")
```
""".format(png=frag_cnt_png)

    print """\n## Base content of sequenced fragments

**ACGT content**: 测序片段中A C G T四种碱基的占比分布。

**CpG content**: 测序片段中CpG二碱基的占比分布, 计算方式为 $\\frac{CpG*2}{(sequence_length)}$

**Observed to expected CpG ratio**: 实际CpG含量与期望CpG含量的比值，计算方式为$Obs/Exp CpG = \\frac{Number of CpG * N}{(Number of C * Number of G)}$ (N = length of sequence)

**CpG island**: 长度200 bp；GC比例50%以上；observed-to-expected CpG ratio 大于0.6。

"""

    frag_base_pdf = report_sub_dir+'/'+os.path.split(frag_base_pdf)[-1]
    frag_base_png = frag_base_pdf.replace('pdf', 'png')

    print "\n(ref:frag-base-content-fig) 测序测到的DNA序列片段的碱基构成。**ACGT content**: 测序片段中A C G T四种碱基的占比分布; **CpG content**: 测序片段中CpG二碱基的占比分布;**Observed to expected CpG ratio**: 实际CpG含量与期望CpG含量的比值。[PDF]({})\n".format(frag_base_pdf)

    print """```{{r frag-base-content-fig, fig.cap="(ref:frag-base-content-fig)"}}
knitr::include_graphics("{png}")
```
""".format(png=frag_base_png)
    if appendL:
        appendList = generateLink(appendL, labelL, 'seq_fragment_attribute', 
                report_sub_dir, '\n* ')
        print """## 测序片段统计表格下载 {#seq-frag-xls}

* {}
""".format(appendList)
#-------------------------------

def output(aDict, label, file_fh, set):
    for key, value in aDict.items():
        print >>file_fh, "{}\t{}\t{}\t{}".format(key, label, value, set)
#---------------------------------------------------------------

def output_base(baseD, frag_base_fh, label):
    #"variable\tvalue\ttype\tset"
    typeD = {'A':'ACGT content (%)', 'T':'ACGT content (%)', 'C':'ACGT content (%)', 
            'G':'ACGT content (%)', 'CpG':'CpG content (%)', 
            'CpG_ratio':'Observed to expected CpG ratio'}
    for key, valueL in baseD.items():
        key = key.replace('_percent', '')
        if key == 'CpG_observe_to_expect_ratio':
            key = "CpG_ratio"
        type = typeD[key]
        for value in valueL:
            print >>frag_base_fh, '{}\t{}\t{}\t{}'.format(key, value, label, type)
#---------------------------------------------------------------
def plot2(frag_len, frag_cnt):
    cmd = ["s-plot lines -m TRUE -a frag_len -V \'c(150)\'", 
            "-x \"Fragment size (bp)\" -y \"Frequency (%) [lower subplot] or raw counts [upper subplot]\" ", 
            "-f", frag_len, '-F \"+facet_wrap(~set, ncol=1, scale=\'free_y\')\" -u 15']
    #print ' '.join(cmd)
    os.system(' '.join(cmd))
    
    cmd = ["s-plot lines -m TRUE -a frag_len -V \'c(150)\' -o TRUE", 
            "-x \"Fragment size (bp)\" -y \"Frequency (%) [lower subplot] or raw counts [upper subplot]\" ", 
            "-f", frag_len, '-F \"+facet_wrap(~set, ncol=1, scale=\'free_y\')\" -u 15']
    #print ' '.join(cmd)
    os.system(' '.join(cmd))

    cmd = ["s-plot lines -m TRUE -a frag_cnt ", 
            "-x \"Fragment occurence time\" -y \"Frequency (%) [lower subplot] or raw counts [upper subplot]\" ", 
            "-f", frag_cnt, '-F \"+facet_wrap(~set, ncol=1, scale=\'free_y\')\" -u 15']
    #print ' '.join(cmd)
    os.system(' '.join(cmd))
#---------------------------------------------------
def plotBase(frag_base):
    #"variable\tvalue\tlabel\tset"
    cmd = ['s-plot boxplot -m TRUE', '-a label -b 45 -W TRUE',
           '-y "Base content"', '-u 25', 
           '-l "\'A\',\'C\',\'T\',\'G\',\'CpG\',\'CpG_ratio\'"', 
           '-f', frag_base, '-p "+facet_wrap(~set,ncol=1,scale=\'free_y\')"']
    print ' '.join(cmd)
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
    frag_len = op+'.fragment_len_frequency.xls'
    frag_cnt  = op+'.fragment_cnt_frequency.xls'
    frag_base = op+'.fragment_base_frequency.xls'

    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    doc_only = options.doc_only
    curation_label = os.path.split(sys.argv[0])[-1].replace('.',  '_')
    if doc_only:
        generateDoc(report_dir, report_sub_dir, frag_len, frag_cnt, frag_base, 
                curation_label, appendL, labelL)
        return 0
    #op_fh = open(op, 'w')
    global debug
    debug = options.debug
    #-----------------------------------
    #i = -1
    #print >>op_fh, "variable\tvalue\tset"
    #for file in fileL:
    #    i+=1
    #    label = labelL[i]
    #    readFile(file, label, op_fh)
    #op_fh.close()

    #plot(op)

    frag_len_fh = open(frag_len, 'w')
    frag_cnt_fh = open(frag_cnt, 'w')
    frag_base_fh = open(frag_base, 'w')

    print >>frag_len_fh, "frag_len\tvariable\tvalue\tset"
    print >>frag_cnt_fh, "frag_cnt\tvariable\tvalue\tset"
    print >>frag_base_fh, "variable\tvalue\tlabel\tset"
    
    i = -1
    
    for file in fileL:
        i+=1
        label = labelL[i]
        lenD, lenFreqD, cntD, cntFreqD, baseD = cntFile(file, test)
        
        output(lenD, label, frag_len_fh, 'Fragment length distrib (count)')
        output(lenFreqD, label, frag_len_fh, 'Fragment length distrib (frequency (%))')

        output(cntD, label, frag_cnt_fh, 'Fragment occurence number distrib (count)')
        output(cntFreqD, label, frag_cnt_fh, 'Fragment occurence number distrib (frequency (%))')
        
        output_base(baseD, frag_base_fh, label)
    frag_len_fh.close()
    frag_cnt_fh.close()
    frag_base_fh.close()

    plot2(frag_len, frag_cnt)
    plotBase(frag_base)

    generateDoc(report_dir, report_sub_dir, frag_len, frag_cnt, frag_base, 
            curation_label, appendL, labelL)

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


