#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, with_statement
'''
Copyright 2013, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Functional description:
    This is used to statistics bismark mapping quality.
'''

import sys
import re
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from tools import *

def fprint(content):
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
    parser.add_option("-i", "--input-files", dest="filein",
        metavar="FILEIN", help="The output of bismark \
bismark_bt2_PE_report.txt. \
Multiple files can be separated with \
<,> or space(< >) or both, like <'file1,file2,file3'>, \
<'file1 file2 file3'> or <'file1' 'file2', 'file3'>." )
    parser.add_option("-l", "--file-label", dest="file_label",
        metavar="FILE_LABEL", help="The labels should have the same \
format and order as filenames given to <-i>.")
    parser.add_option("-S", "--seq-type", dest="seq_type",
        default='PE', help="PE or SE")
    parser.add_option("-o", "--output-prefix", dest="op",
        help="Specify output file prefix.")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='2_mapping_quality', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '2_mapping_quality'.")
    parser.add_option("-a", "--appendFile", dest="append",
        default='', help="A list of files to be appended. Unused.")
    parser.add_option("-d", "--doc-only", dest="doc_only",
        default=False, action="store_true", help="Specify to only generate doc.")
    parser.add_option("-t", "--test-line", dest="test",
        type="int", default=0, help="Give a positive number (n) to test the program using on the first n reads. Default 0 meaning using all reads in file.")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readBisulfite(file, seq_type):
    mapD = {'Sequence pairs analysed in total': 'Total read pairs', 
            'Number of paired-end alignments with a unique best hit': 
            'Unique-mapped read pairs', 
            'Sequence pairs with no alignments under any condition': 
            'Unmapped read pairs', 
            'Sequence pairs did not map uniquely':'Multi-mapped read pairs', 
            "Total number of C's analysed": "Total analysed C", 
            "Total methylated C's in CpG context": "mCpG", 
            "Total methylated C's in CHG context": "mCHG", 
            "Total methylated C's in CHH context": "mCHH", 
            "Total methylated C's in Unknown context": "mCother", 
            "Total unmethylated C's in CpG context": "umCpG", 
            "Total unmethylated C's in CHG context": "umCHG", 
            "Total unmethylated C's in CHH context": "umCHH", 
            "Total unmethylated C's in Unknown context": "umCother"
            }
    resultDict = {}
    if seq_type != 'PE':
        print >>sys.stderr, "Not supported for single end sequencing"
        sys.exit(1)
    #-----------------------------------
    begin = 1
    for line in open(file):
        #if begin and not line.startswith("Final Alignment report"):
        #    continue
        #else:
        #    begin = 0
        if line.find(':') != -1:
            key, value = line.split(':')[:2]
            if key in mapD:
                newKey = mapD[key]
                value = int(value.strip().split()[0])
                resultDict[newKey] = value
            if key == "Sequence pairs which were discarded because genomic sequence could not be extracted":
                value = int(value.strip().split()[0])
                resultDict['Unmapped read pairs'] += value
    assert len(resultDict) == len(mapD), resultDict
    return resultDict
#-----------------------------------------------------------------------------
def plot(read_melt_file, c_melt_file, nameL):
    '''
    variable value sample set
    '''
    x_level = ["'"+i+"'" for i in nameL]
    x_level = '"'+','.join(x_level)+'"'
    cmd = ["s-plot barPlot -m TRUE -a sample -R 90 -B set -O 1 -w 20 -u 20 -f ", 
            read_melt_file, '-P top -k free_y -L', x_level, 
            ' -y  \'Reads count or relative percent\' -x \'Samples\' ']
    #print ' '.join(cmd)
    os.system(' '.join(cmd))
    cmd = ["s-plot barPlot -m TRUE -a sample -R 90 -B set -O 1 -w 20 -u 20 -f ", 
            c_melt_file, ' -k free_y -L', x_level, 
            "-l \"'mCpG','mCHG','mCHH','mCother','umCpG','umCHG','umCHH','umCother'\"", 
            "-c TRUE -C 'colorRampPalette(c(\"red\",\"orange\",\"dark green\"),alpha=0.5)(8)'", 
            ' -y  \'Raw count or relative percent\' -x \'Samples\' ']
    #print >>sys.stderr, ' '.join(cmd)
    os.system(' '.join(cmd))
#--------------------------------------------------------

def generateDoc(report_dir, report_sub_dir, table_file, read_melt_file,  
        c_melt_file, curation_label, labelL, melt=1):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)
    if melt:
        read_pdf = read_melt_file+'.stackBars.pdf'
        c_pdf = c_melt_file+'.stackBars.pdf'
        copypdf(dest_dir, read_pdf, c_pdf)
    else:
        per_pdf = totalTable_per+'.heatmapS.pdf'
        cnt_pdf = totalTable_cnt+'.heatmapS.pdf'
        copypdf(dest_dir, per_pdf, cnt_pdf)
    copy(dest_dir, table_file, read_melt_file, c_melt_file)

    print "\n## Tabulated statistics of reads mapping\n"

    knitr_read_txt(report_dir,  curation_label)
    
    print """每个样品测序序列比对总结见 (Table \@ref(tab:read-map-sum-cnt-table) and Figure \@ref(fig:read-map-sum-fig))。

```{r describe-tab-map-status}
text = "Label;Explanation
Total read pairs;预处理后用于比对的reads数 
Unique-mapped read pairs;在基因组上有唯一比对位置的reads数
Multi-mapped read pairs;在基因组上有多个比对位置的reads数
Total analysed C;测序覆盖的胞嘧啶的数目
mCpG;甲基化的CpG的数目
mCHG;CHG序列中甲基化的C的数目，H代表非G的脱氧核糖核苷酸
mCHH;CHH序列中甲基化的C的数目，H代表非G的脱氧核糖核苷酸
mCother;其它序列中甲基化的C的数目
umCpG;未甲基化的CpG的数目
umCHG;CHG序列中未甲基化的C的数目，H代表非G的脱氧核糖核苷酸
umCHH;CHH序列中未甲基化的C的数目，H代表非G的脱氧核糖核苷酸
umCother;其它序列中未甲基化的C的数目"

describe_tab <- read.table(text=text, sep=";", header=1, quote="")
knitr::kable(describe_tab)
```

"""

    print "Reads比对数据的统计可以帮助判断测序的质量、准确度、有无偏好和异常等。如果过滤后用于下游分析的reads数偏少，则需要慎重考虑。对于未比对回去的reads需要考虑未比对回去的原因，区别对待。\n"

    print "原始数据或者PDF格式的文件可以点击XLS或PDF下载。\n"
    
    totalTable_cntNew = report_sub_dir+'/'+os.path.split(table_file)[-1]
    print "(ref:read-map-sum-cnt-table) Summary raw counts of mapped reads. [XLS]({})\n".format(totalTable_cntNew)
    
    print '''```{{r read-map-sum-cnt-table}}
map_table <- read.table("{totalTable_cntNew}", sep="\\t", header=T, row.names=1, quote="", comment="")
knitr::kable(map_table, booktabs=T, caption="(ref:read-map-sum-cnt-table)")
```
'''.format(totalTable_cntNew=totalTable_cntNew)

    #totalTable_perNew = report_sub_dir+'/'+os.path.split(totalTable_per)[-1]
    #print "(ref:read-map-sum-per-table) Summary percentage of mapped reads relative to total reads. [XLS]({})\n".format(totalTable_perNew)
    
    #print '''```{{r read-map-sum-per-table}}
#map_table <- read.table("{totalTable_perNew}", sep="\\t", header=T, row.names=1, quote="", comment="")
#knitr::kable(map_table, booktabs=T, caption="(ref:read-map-sum-per-table)")
#```
#'''.format(totalTable_perNew=totalTable_perNew)
    
    if melt:
        #read_pdf = read_melt_file+'.stackBars.pdf'
        #c_pdf = c_melt_file+'.stackBars.pdf'
        pdf = report_sub_dir+'/'+os.path.split(read_pdf)[-1]
        xls = report_sub_dir+'/'+os.path.split(read_melt_file)[-1]
        print "(ref:read-map-sum-fig) Summary of reads mapping status. \
1 Million = 10^6^. [PDF]({}) [XLS]({})\n".format(pdf, xls)

        png = pdf.replace('pdf', 'png')
        print '''```{{r read-map-sum-fig, fig.cap="(ref:read-map-sum-fig)"}}
knitr::include_graphics("{png}")
```
'''.format(png=png)

        pdf = report_sub_dir+'/'+os.path.split(c_pdf)[-1]
        xls = report_sub_dir+'/'+os.path.split(c_melt_file)[-1]
        print "(ref:c-sum-fig) Summary of C status. \
[PDF]({}) [XLS]({})\n".format(pdf, xls)

        png = pdf.replace('pdf', 'png')
        print '''```{{r c-sum-fig, fig.cap="(ref:c-sum-fig)"}}
knitr::include_graphics("{png}")
```
'''.format(png=png)
    else:
        cnt_pdf = report_sub_dir+'/'+os.path.split(cnt_pdf)[-1]
        cnt_png = cnt_pdf.replace('pdf', 'png')
        per_pdf = report_sub_dir+'/'+os.path.split(per_pdf)[-1]
        per_png = per_pdf.replace('pdf', 'png')

        print "(ref:read-map-sum-fig) Summary of reads mapping status. 1 Million = 10^6^. \
[PDF_cnt]({}) [PDF_percent]({})\n".format(cnt_pdf, per_pdf)

        print '''```{{r read-map-sum-fig, out.width="49%", fig.cap="(ref:read-map-sum-fig)"}}
knitr::include_graphics(c("{cnt_png}", "{per_png}"))
```
'''.format(cnt_png=cnt_png, per_png=per_png)

#--------------------------------
#--------------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    fileL = re.split('[, ]*', options.filein.strip())
    labelL = re.split('[, ]*', options.file_label.strip())
    lenfileL = len(fileL)
    lenlabeL = len(labelL)
    assert lenfileL == lenlabeL
    debug = options.debug
    seq_type = options.seq_type
    table_file = options.op + '.table.xls'
    read_melt_file = options.op + '.read_melt.xls'
    c_melt_file = options.op + '.c_melt.xls'
    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    doc_only = options.doc_only
    curation_label = os.path.split(sys.argv[0])[-1].replace('.',  '_')
    if doc_only:
        generateDoc(report_dir, report_sub_dir, table_file, read_melt_file,  
                c_melt_file, curation_label, labelL)
        return 0
    #-----------------------------------
    table_file_fh = open(table_file, 'w')
    read_melt_file_fh = open(read_melt_file, 'w')
    c_melt_file_fh = open(c_melt_file, 'w')
    keyL = ['Total read pairs', 'Unique-mapped read pairs', 
            'Unmapped read pairs', 
            'Multi-mapped read pairs', "Total analysed C", 
            "mCpG", "mCHG", "mCHH", "mCother", "umCpG", 
            "umCHG", "umCHH", "umCother"]

    setD = {"Reads mapping status":['Unique-mapped read pairs', 
            'Unmapped read pairs', 'Multi-mapped read pairs'], 
            "Different types of Cs":    
            ["mCpG", "mCHG", "mCHH", "mCother", "umCpG", 
            "umCHG", "umCHH", "umCother"]
        }

    print >>table_file_fh, "Sample\t{}".format('\t'.join(keyL))
    print >>read_melt_file_fh, "variable\tvalue\tsample\tset"
    print >>c_melt_file_fh, "variable\tvalue\tsample\tset" 

    for file, label in zip(fileL, labelL):
        resultD = readBisulfite(file, seq_type)
        valueL = [str(resultD[key]) for key in keyL]
        print >>table_file_fh, "{}\t{}".format(label, '\t'.join(valueL))

        key = "Reads mapping status"
        variableL = setD[key]
        total_read = resultD['Total read pairs']/100
        for variable in variableL:
            tmp_151 = [variable, str(resultD[variable]), label, key+' (raw count)']
            print >>read_melt_file_fh, "\t".join(tmp_151)
            tmp_151 = [variable, str(resultD[variable]/total_read), label, 
                key+' (relative percent)']
            print >>read_melt_file_fh, "\t".join(tmp_151)

        key = "Different types of Cs"
        variableL = setD[key]
        total_c = resultD["Total analysed C"]/100
        for variable in variableL:
            tmp_151 = [variable, str(resultD[variable]), label, key+' (raw count)']
            print >>c_melt_file_fh, "\t".join(tmp_151)
            tmp_151 = [variable, str(resultD[variable]/total_c), label, key+' (relative percent)']
            print >>c_melt_file_fh, "\t".join(tmp_151)

    table_file_fh.close()
    read_melt_file_fh.close()
    c_melt_file_fh.close()
    plot(read_melt_file, c_melt_file, labelL)

    generateDoc(report_dir, report_sub_dir, table_file, read_melt_file,  
            c_melt_file, curation_label, labelL)
#--------------main-----------------------------

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


