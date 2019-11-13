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
    This is designed to summarize reads distribution output by `read_distribution.py`.
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
        metavar="FILEIN", help="`,` or ` ` separated a list of files.")
    parser.add_option("-l", "--labels", dest="label",
        metavar="LABEL", help="`,` or ` ` separated a list of labels to label each file. It must have same order as files.")
    parser.add_option("-c", "--chromosome-size", dest="chr_size",
        help="Chromsome size file. A tab separated file with first column as chr_name and second column as chr_size.")
    parser.add_option("-o", "--output-file", dest="output",
        help="Name for output file.")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='2_mapping_quality', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '2_mapping_quality'.")
    parser.add_option("-d", "--doc-only", dest="doc_only",
        default=False, action="store_true", help="Specify to only generate doc.")
    parser.add_option("-n", "--number", dest="number", type="int", 
        default=40, help="Set the maximum allowed samples for barplot. Default 40.\
 If more than this number of samples are given, heatmap will be used instead. ")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readFile(file, genome):
    fh = open(file)
    total_reads = int(fh.readline().rsplit()[-1])
    total_tags = int(fh.readline().rsplit()[-1])
    total_assigned_tags = int(fh.readline().rsplit()[-1])
    split_line = fh.readline()
    assert split_line.startswith('===='), split_line
    header = fh.readline()
    countD = {}
    percentD = {}
    normalizeCountD = {}
    baseD = {}
    
    while 1:
        line = fh.readline()
        if not line or line.startswith('==='):
            break
        #print >>sys.stderr, line, 
        type, base, count, normC = line.split()
        type = type.replace('_Exons', '')
        countD[type] = int(count)
        normalizeCountD[type] = float(normC)
        baseD[type] = int(base)
    #--------------------------
    fh.close()
    classT = ['TSS_up', 'TES_down']
    posT   = ['10kb', '5kb', '1kb']
    
    for c in classT:
        for i in range(2):
            pos = posT[i]
            pos_before = posT[i+1]
            name = c+'_'+pos
            name_before = c+'_'+pos_before
            countD[name] = countD[name]-countD[name_before]
            baseD[name] = baseD[name]-baseD[name_before]
            normalizeCountD[name] = countD[name]/baseD[name]* 1000

    countD['Intergenic'] = total_tags - total_assigned_tags
    normalizeCountD['Intergenic'] = countD['Intergenic']*1000/(genome-sum(baseD.values()))
    sumall = sum(countD.values())
    assert sumall == total_tags
    for key, value in countD.items():
        percentD[key] = 100.0 * value/sumall

    #--------------------
    normalizepercentD = {}
    sum_normalize = sum(normalizeCountD.values())
    for key, value in normalizeCountD.items():
        normalizepercentD[key] = 100.0 * value/sum_normalize


    return countD, percentD, normalizeCountD, normalizepercentD
#-, ------------------------

def plot(op, nameL):
    x_level = ["'"+i+"'" for i in nameL]
    x_level = '"'+','.join(x_level)+'"'
    cmd = "s-plot barPlot -m TRUE -a Sample -R 45 -B type -O 1 -w 20 -u 30 -f " \
            + op + ' -x \'Samples\' -k free_y' + ' -L ' + x_level
    os.system(cmd)
#------------------------------------

def plot_heatmap(op, typeL, sample_cnt):
    '''
    typeL = ['Read_cnt','Read_percent','Read_cnt_normalize','Read_percent_normalize']
    '''
    height = sample_cnt // 2 + 5
    if height < 10:
        height = 10
    for type in typeL:
        file = op+'.'+type+'.xls'
        cmd = ["s-plot heatmapS -a TRUE -A 45 -b TRUE -R TRUE", 
                "-x white -y blue -u 17 -F 12 -v", str(height), 
                "-f ", file, "-I", type, "-l top -T 2", 
                ]
        os.system(' '.join(cmd))
#----------------------------------------------

def generateDoc(report_dir, report_sub_dir, op, curation_label, heatmap):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)

    print "\n## Reads distribution along genome landmark regions\n"
    
    curation_label = "Reads_distribution_along_genome_landmark_regions"
    knitr_read_txt(report_dir,  curation_label)

    print """基因组标志性区域的reads分布情况统计。

转录组测序的是转录的RNA，通常情况下是CDS, 3'UTR, 5'UTR reads最多；Intron区的reads表明存在`Intron retention`类型的可变剪切；Intergenic区的reads表明存在未注释的基因。

1. **Read_cnt** 表示落在对应基因组区间的reads数；
2. **Read_cnt_normalize** 表示每个基因组区间平均每KB区域的reads数. [从分布随机性来看，越长的基因组区域就可能统计出越多的READs分布, 因此需要计算单位长度READs的数目以便比较READs在某个基因组区域的富集情况。]；
3. **Read_percent** 表示落在对应基因组区间的reads的相对比例；计算这个是为了方便样品之间的比较，虽然测序深度不同，可以比较样品间落在特定区域的reads的比例的差异。
4. **Read_percent_normalize** 表示每个基因组区间平均每KB区域的reads的相对比例。
  
区域解释：

* TSS_up_1kb: 转录起始位点上游 1 kb; TSS_up_5kb：转录起始位点上游 5 kb (不包含前面提到的1 kb)。
* TES_down_1kb: 转录终止位点下游 1 kb; TES_down_5kb：转录终止位点下游 5 kb (不包含前面提到的1 kb)。

"""
    if heatmap:
        typeL = ['Read_cnt','Read_percent','Read_cnt_normalize',
                'Read_percent_normalize']
        tableL = [op+'.'+type+'.xls' for type in typeL]
        pdfL = [table+'.heatmapS.pdf' for table in tableL]
        pngL = [table+'.heatmapS.png' for table in tableL]
        copy(dest_dir, *tableL)
        copypdf(dest_dir, *pdfL)

        
        #pdf_link = "[label_pdf](pdf), [label_pdf](pdf)"
        pdf_link = generateLink(pdfL, typeL, 'pdf', report_sub_dir)
        xls_link = generateLink(tableL, typeL, 'table', report_sub_dir)

        print "(ref:reads-distrib-sum-fig) Reads distribution along genome landmark regions. **Read_cnt** 表示落在对应基因组区间的reads数；**Read_cnt_normalize** 表示每个基因组区间平均每KB区域的reads数；**Read_percent**表示落在对应基因组区间的reads的相对比例；**Read_percent_normalize** 表示每个基因组区间平均每KB区域的reads的相对比例。 {} {}\n".format(pdf_link, xls_link)
        png = grenerateQuotedLists(pngL, report_sub_dir)
        print """```{{r reads-distrib-sum-fig, out.width="49%", fig.cap="(ref:reads-distrib-sum-fig)"}}
knitr::include_graphics(c({png}))
```
""".format(png=png)
    else:
        totalTable = op + '.xls'
        pdf = totalTable+'.stackBars.pdf'
        copy(dest_dir, totalTable)
        copypdf(dest_dir, pdf)
        pdf = report_sub_dir+'/'+os.path.split(pdf)[-1]
        newTable = report_sub_dir+'/'+os.path.split(totalTable)[-1]
        png = pdf.replace('pdf', 'png')

        print "(ref:reads-distrib-sum-fig) Reads distribution along genome landmark regions. **Read_cnt** 表示落在对应基因组区间的reads数；**Read_cnt_normalize** 表示每个基因组区间平均每KB区域的reads数；**Read_percent**表示落在对应基因组区间的reads的相对比例；**Read_percent_normalize** 表示每个基因组区间平均每KB区域的reads的相对比例。 [PDF]({}) [XLS]({})\n".format(pdf, newTable)

        print """```{{r reads-distrib-sum-fig, fig.cap="(ref:reads-distrib-sum-fig)"}}
knitr::include_graphics("{png}")
```
""".format(png=png)
#-------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    fileL = re.split(r'[, ]*', file.strip())
    label = options.label
    labelL = re.split(r'[, ]*', label.strip())
    chr_size = options.chr_size
    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    doc_only = options.doc_only
    op = options.output
    sample_cnt = len(fileL)
    allowed_sample = options.number
    if sample_cnt > allowed_sample:
        heatmap = 1
    else:
        heatmap = 0
    curation_label = os.path.split(sys.argv[0])[-1].replace('.', '_')
    if doc_only:
        generateDoc(report_dir, report_sub_dir, op, curation_label, heatmap)
        return 0
    genome = sum([int(line.split()[1]) for line in open(chr_size)])
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    aDict = {}
    if heatmap:
        headerL = ['TSS_up_10kb', 'TSS_up_5kb', 'TSS_up_1kb', '5\'UTR', 'CDS', 
                'Introns', '3\'UTR', 'TES_down_1kb', 'TES_down_5kb', 'TES_down_10kb', 
                'Intergenic']
        fileD_fh = {"Read_cnt": open(op+'.Read_cnt.xls', 'w'), 
            "Read_percent": open(op+'.Read_percent.xls', 'w'), 
            "Read_cnt_normalize": open(op+'.Read_cnt_normalize.xls', 'w'), 
            "Read_percent_normalize": open(op+'.Read_percent_normalize.xls', 'w')
            }
        for fh_label, fh in fileD_fh.items():
            print >>fh, "sample\t{}".format('\t'.join(headerL))
    else:
        op_fh = open(op+'.xls', 'w')
        print >>op_fh, "Sample\tvariable\tvalue\ttype"
    i = -1
    for file in fileL:
        i += 1
        label = labelL[i]
        countD, percentD, normalizeCountD, normalizepercentD = readFile(file, genome)
        tmpD = {'Read_cnt': countD, 'Read_percent': percentD, 
                'Read_cnt_normalize': normalizeCountD, 
                'Read_percent_normalize': normalizepercentD}
        if heatmap:
            for fh_label, fh in fileD_fh.items():
                dataD = tmpD[fh_label]
                tmpLineL = [str(dataD[region]) for region in headerL]
                print >>fh, "{}\t{}".format(label, '\t'.join(tmpLineL))
        else:
            keyL = countD.keys()
            keyL.sort()
            for type, typeD in tmpD.items():
                for key in keyL:
                    print >>op_fh, '{}\t{}\t{}\t{}'.format(
                     label, key, typeD[key], type)
                ###--------multi-process------------:------
    #-----------------------------------------
    if heatmap:
        for fh in fileD_fh.values():
            fh.close()
        typeL = ['Read_cnt','Read_percent','Read_cnt_normalize','Read_percent_normalize']
        plot_heatmap(op, typeL, sample_cnt)
    else:
        op_fh.close()
        plot(op+'.xls', labelL)
    generateDoc(report_dir, report_sub_dir, op, curation_label, heatmap)

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


