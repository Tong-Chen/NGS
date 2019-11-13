#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is used to generate FASTQC report.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
from tools import *
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

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
        metavar="FILEIN", help="Space separated sample names")
    parser.add_option("-s", "--seq-type", dest="seqTypeL",
        metavar="seqTypeL", default="makefile.am.template", 
        help="Default makefile.am.template. \
Only the line containing `seq_type=` will be used.")
    parser.add_option("-d", "--dir", dest="dir",
        default="report", help="The dir for document and result files.\
default <report>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_false", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def qualityD_plot(qualityD, main_dir, file_relative):
    '''
    qualityD = {"Sample" : {"GC": "PASS", "base":"WARN"}}
    '''
    #file_relative = tar_dir + '/fastqc_summary.txt'
    file = main_dir + '/' + file_relative
    fh = open(file, 'w')
    print >>fh, "Sample\tGC_quality\tBase_quality"
    for sample, subD in qualityD.items():
        print >>fh, "{}\t{}\t{}".format(sample, subD["Per sequence GC content"], subD["Per base sequence quality"])
    fh.close()
    return file_relative
#-----------------------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    curation_label = os.path.split(sys.argv[0])[-1].replace('.', '_')
    main_dir = options.dir
    sampleL = options.filein.split()
    seqTypeL = options.seqTypeL
    seqTypeD = {}
    for line in open(seqTypeL):
        if line.find("_seq_type=") != -1:
            samp, seqType = line.strip().split('_seq_type=')
            seqTypeD[samp] = seqType
    #----------------------------------------------
    tardir = "1_sequencing_quality_check"
    os.system("mkdir -p "+main_dir+'/'+tardir)
    
    file_relative = tardir + '/fastqc_summary.txt'
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    print """## 测序碱基质量、接头状态、GC含量和均衡度
"""
    curation_label = "Fastqc_result_summary"
    knitr_read_txt(main_dir, curation_label)
    
    print """

测序质量概览如 Figure \@ref(fig:fastqc-summary-scatter-plot) 所示。`GC_quality`检测测序数据的GC含量是否正常。`Base_quality`检测测序数据的碱基偏好性。如果是Bisulfite测序，请忽略这个结果。

```{{r fastqc-summary-scatter-plot, cache.extra=mtime("{file}"), fig.height=9, fig.width=9, fig.cap="Summary of sequencing quality. GCquality (Per sequence GC content); Basequality (Per base sequence quality)"}}
library(ggplot2)
library(ggbeeswarm)
data <- read.table("{file}", header=T, sep="\t", quote="")
ggplot(data, aes(x=GC_quality, y=Base_quality)) + geom_quasirandom(aes(color=GC_quality, shape=Base_quality), groupOnX=FALSE, alpha=0.7) + theme_bw() + theme(panel.grid.minor = element_blank()) + geom_text(aes(label=Sample), position=position_quasirandom(groupOnX=F), size=2)
```

""".format(file=file_relative)

    print """
Table: (\#tab:seq-quality-explanatioan-ch) 测序质量评估结果解读方法

-----------------------------------------------------------------------------------
评估内容                   结果解释 (图例中会标记对应评估内容为PASS、WARN和FAIL, 具体处理方式详见下面中英文解释)
-------------------------  --------------------------------------------------------------------------------
Per base quality           测序reads从5'到3'的碱基的质量值 (Q)。该值越大越代表对应碱基测序准确度越高。假设p为一个碱基测序错误的概率，则Q=-10 * log10(p). 质量值为10时，对应碱基出错的概率为10%；质量值为20时，对应碱基出错的概率为1%。通常来讲，3'端的碱基质量会低于5'端；另外5'端最初几个碱基也会出现较大的质量值波动。我们在后期处理时，会去除低质量的碱基以保证分析结果的准确性。

Adaptor content            判断测序reads中是否残留接头序列。存在接头序列和不存在接头序列都是合理的，取决于测序数据下机后是否进行了接头去除和去除的是否完整。若在分析时检测到接头序列存在，我们会首先去除接头，然后进行后续分析，以保证分析结果的准确性。

Per sequence GC content    测序reads的GC含量。正常的测序reads的GC含量符合正态分布模式 (形如图中蓝色的倒钟形线)。若在平滑的曲线上存在一个尖峰表示测序样品存在特定的序列污染如混入了引物二聚体。若GC含量分布曲线比较平坦则代表可能存在不同物种的序列污染。当这一指标异常时，可能导致后期的序列比对或拼接存在问题，需要引起注意。

Per base sequence content  测序reads的碱基偏好性。正常的测序结果中一个序列不同的碱基没有偏好性，图中的线应平行。Bisulfite测序中存在甲基化的C到T的转变，会导致这一评估结果异常。我们人工核验无误后，可以忽略软件对这一检测结果的评价。
-----------------------------------------------------------------------------------


Table: (\#tab:seq-quality-explanatioan-en) Explanation of each picture.

-----------------------------------------------------------------------------------
Analysis                   Explanation
-------------------------  --------------------------------------------------------------------------------
Per base quality           The most common reason for warnings and failures in this module is a general degradation of quality over the duration of long runs. In general sequencing chemistry degrades with increasing read length and for long runs you may find that the general quality of the run falls to a level where a warning or error is triggered.

Per sequence GC content    Warnings in this module usually indicate a problem with the library. Sharp peaks on an otherwise smooth distribution are normally the result of a specific contaminant (adapter dimers for example),  which may well be picked up by the overrepresented sequences module. Broader peaks may represent contamination with a different species.

Adaptor content            Any library where a reasonable proportion of the insert sizes are shorter than the read length will trigger this module. This doesn't indicate a problem as such - just that the sequences will need to be adapter trimmed before proceeding with any downstream analysis.

Per base sequence content  In a random library you would expect that there would be little to no difference between the different bases of a sequence run,  so the lines in this plot should run parallel with each other. The relative amount of each base should reflect the overall amount of these bases in your genome,  but in any case they should not be hugely imbalanced from each other.
-----------------------------------------------------------------------------------


"""
    label = 0
    qualityD = {}
    for sample in sampleL:
        label += 1
        if seqTypeD[samp] == 'SE':
            summary_file = sample+"_fastqc/summary.txt"
            summaryD = dict([(line.split('\t')[1],line.split('\t')[0]) \
                for line in open(summary_file)])
            assert sample not in qualityD, "Duplicate sample "+sample
            qualityD[sample] = {}
            qualityD[sample]["Per sequence GC content"] = summaryD["Per sequence GC content"]
            qualityD[sample]["Per base sequence quality"] = summaryD["Per base sequence quality"]
            os.system("/bin/cp -u --parent %s_fastqc/Images/* %s/%s/" \
                % (sample, main_dir, tardir))
            print """### Sequencing quality estimation for sample *{sample}*

(ref:seq-quality-{label}) These 4 pictures show (top-left to bottom-right): Per base quality *({summaryD[Per base sequence quality]})*, Adaptor content *({summaryD[Adapter Content]})*, Per sequence GC content *({summaryD[Per sequence GC content]})*, Per base sequence content *({summaryD[Per base sequence content]})*.

```{{r seq-quality-{label}, fig.cap="(ref:seq-quality-{label})", out.width="49%"}}
figs_{label} = paste0("{dir}", c("per_base_quality.png", "adapter_content.png", "per_sequence_gc_content.png", "per_base_sequence_content.png"))
knitr::include_graphics(figs_{label})
```
""".format(label=label, summaryD=summaryD, dir=tardir+'/'+sample+'_fastqc/Images/', sample=sample)

        elif seqTypeD[samp] == 'PE':
            summary_file = sample+"_1_fastqc/summary.txt"
            summaryD1 = dict([(line.split('\t')[1],line.split('\t')[0]) \
                for line in open(summary_file)])
            os.system("/bin/cp -u --parent %s_1_fastqc/Images/* %s/%s/" \
                % (sample, main_dir, tardir))
            summary_file = sample+"_2_fastqc/summary.txt"
            summaryD2 = dict([(line.split('\t')[1],line.split('\t')[0]) \
                for line in open(summary_file)])
            os.system("/bin/cp -u --parent %s_2_fastqc/Images/* %s/%s/" \
                % (sample, main_dir, tardir))
            
            assert sample+'_1' not in qualityD, "Duplicate sample "+sample
            qualityD[sample+"_1"] = {}
            qualityD[sample+"_1"]["Per sequence GC content"] = summaryD1["Per sequence GC content"]
            qualityD[sample+"_1"]["Per base sequence quality"] = summaryD1["Per base sequence quality"]
            qualityD[sample+"_2"] = {}
            qualityD[sample+"_2"]["Per sequence GC content"] = summaryD2["Per sequence GC content"]
            qualityD[sample+"_2"]["Per base sequence quality"] = summaryD2["Per base sequence quality"]
            print """### Sequencing quality estimation for sample *{sample}*

(ref:seq-quality-{label}) These 8 pictures show (top-left to bottom-right): Per base quality for left reads *({summaryD1[Per base sequence quality]})* and right reads *({summaryD2[Per base sequence quality]})*, Adaptor content for left reads *({summaryD1[Adapter Content]})* and right reads *({summaryD2[Adapter Content]})*, Per sequence GC content for left reads *({summaryD1[Per sequence GC content]})* and right reads *({summaryD2[Per sequence GC content]})*, Per base sequence content for left reads *({summaryD1[Per base sequence content]})* and right reads *({summaryD2[Per base sequence content]})*.

```{{r seq-quality-{label}, fig.cap="(ref:seq-quality-{label})", out.width="49%"}}
figs_{label} = c(paste0(c("{dir1}", "{dir2}"), "per_base_quality.png"), 
                 paste0(c("{dir1}", "{dir2}"), "adapter_content.png"), 
                 paste0(c("{dir1}", "{dir2}"), "per_sequence_gc_content.png"), 
                 paste0(c("{dir1}", "{dir2}"), "per_base_sequence_content.png")) 
knitr::include_graphics(figs_{label})
```
""".format(label=label, summaryD1=summaryD1, summaryD2=summaryD2, 
        dir1=tardir+'/'+sample+'_1_fastqc/Images/', 
        dir2=tardir+'/'+sample+'_2_fastqc/Images/', sample=sample)

    qualityD_plot(qualityD, main_dir, file_relative)
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


