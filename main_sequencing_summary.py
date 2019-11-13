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
    Read the Fastqc output to get # of sequencing reads, GC content.
    Read the original fastq file to get # of sequenced bases.

---SampleFile-------------------
Samp    conditions
T0_1    T0
T0_2    T0
T0_3    T0
T2_1    T2
T2_2    T2
T2_3    T2



'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from getTopDEgenes import transferListToMultiLTable
#from multiprocessing.dummy import Pool as ThreadPool
from tools import copypdf, copy, knitr_read_txt, transferListToMultiLTable

debug = 0

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
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="Space sepated sample list")
    parser.add_option("-m", "--makefile-am-template",
        dest="makefile_am", default="makefile.am.template", 
        metavar="makefile.am.template", 
        help="Default <makefile.am.template>. \
Only the seq_type information will be used.")
    parser.add_option("-s", "--sampleFile", dest="sampleFile",
        metavar="sampleFile", 
        default="sampleFile", help="Default <sampleFile> with \
format speified above.")
    parser.add_option("-d", "--dir", dest="dir",
        default="report", help="The dir for document and result files.\
default <report>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=False, help="Display running status")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def parseFile(file, subD):
    for line in open(file):
        if line.startswith('>>END_MODULE'):
            break
        key, value = line.strip().split('\t')
        subD[key] = value
#-------------------------------

def readSample(sampleFile):
    '''
    Samp    conditions
    T0_1    T0
    T0_2    T0
    T0_3    T0
    T2_1    T2
    T2_2    T2
    T2_3    T2

    sampD = {'T0_1':"T0", 'T0_2':"T0", 'T0_3':"T0", }

    condD = {"T0": ["T0_1", "T0_2", "T0_3"]}

    '''
    sampD = {}
    condD = {}
    header = 1
    for line in open(sampleFile):
        if header:
            header -= 1
            continue
        #---------------------------
        samp, cond = line.split()[:2]
        assert samp not in sampD, "Duplicated %s" % samp

        sampD[samp] = cond

        if cond not in condD:
            condD[cond] = []
        condD[cond].append(samp)

    return sampD, condD
#--------------------------------

def readBases(read_base_countL, labelMapD):
    baseCntD_T0_1 = {}
    baseCntD_T0_1_1 = {}
    for bcF in read_base_countL:
        for line in open(bcF):
            T0_1_1, type, count = line.strip().split()
            if type == "bases_count":
                assert T0_1_1 not in baseCntD_T0_1_1, "Duplicate "+T0_1_1
                baseCntD_T0_1_1[T0_1_1] = count
                T0_1 = labelMapD[T0_1_1]
                if T0_1 not in baseCntD_T0_1:
                    baseCntD_T0_1[T0_1] = int(count)
                else:
                    baseCntD_T0_1[T0_1] += int(count)
        #-------------------------------------
    #-------------------------------------
    return baseCntD_T0_1, baseCntD_T0_1_1
#---------readBases-----------------



def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    labelL      = options.filein.split(' ')
    makefile_am = options.makefile_am
    sampleFile  = options.sampleFile
    main_dir = options.dir
    sampD, condD = readSample(sampleFile)
    curation_label = os.path.split(sys.argv[0])[-1].replace('.', '_')
    seqTypeD = {}
    for line in open(makefile_am):
        if line.startswith("prefix="):
            prefix = line.strip().split('=')[1]
        elif line.find('seq_type=') != -1:
            sample, seq_type = \
                line.strip().replace('_seq_type','').split('=')
            seqTypeD[sample] = seq_type
    newfileL = []
    read_base_countL   = []
    new_labelL = []
    len_labelL = len(labelL)
    labelMapD = {}

    for i in range(len_labelL):
        label = labelL[i]
        seq_type = seqTypeD[label]
        if seq_type == "PE":
            left = label+'_1'
            right = label+'_2'
            new_labelL.extend([left, right])
            newfileL.extend([left+'_fastqc/fastqc_data.txt', 
                right+'_fastqc/fastqc_data.txt'])
            read_base_countL.extend([left+".statistics_fastq_reads_bases", 
                right+".statistics_fastq_reads_bases"])
            labelMapD[left]  = label
            labelMapD[right] = label
        elif seq_type == "SE":
            new_labelL.append(label)
            newfileL.append(label+'_fastqc/fastqc_data.txt')
            read_base_countL.extend([label+".statistics_fastq_reads_bases"])
            labelMapD[label] = label

    fileL = newfileL
    labelL = new_labelL
    
    #print read_base_countL
    baseCntD_T0_1, baseCntD_T0_1_1 = readBases(read_base_countL,labelMapD)
    #print >>sys.stderr, baseCntD_T0_1
    #print >>sys.stderr, baseCntD_T0_1_1
    
    len_labelL = len(labelL)
    #print fileL
    #print labelL
    verbose = options.verbose
    global dbug
    debug = options.debug
    #-----------------------------------
    print "## 测序质量总结 {#sub-sequence-summary}\n"
    
    curation_label = "Count_sequenced_reads_bases"
    knitr_read_txt(main_dir, curation_label)

    tableL = []
    headerL = ["Sample", "Total reads", "Total bases", 
            "Sequence length (nt)", 
            "GC content (%)", "Encoding"]
    tableL.append(headerL)
    #Do not change keyL
    keyL = ["Total Sequences", "Total bases", "Sequence length", "%GC", "Encoding"]
    aDict = {}
    readCntD = {}
    baseCntD = {}

    T0_L = []
    T0_1_L = []

    sampleRepD = {}
    for i in range(len_labelL):
        file = fileL[i]
        label = labelL[i]
        aDict[label] = {}
        parseFile(file, aDict[label])
        aDict[label]["Total bases"] = baseCntD_T0_1_1[label]
        sampleL = [label]
        sampleL.extend([aDict[label].get(key) for key in keyL])
        tableL.append(sampleL)
        T0_1_1 = label
        T0_1 = labelMapD[T0_1_1]
        #print >>sys.stderr, T0_1_1, T0_1
        T0 = sampD[T0_1]
        #print >>sys.stderr, T0_1_1, T0_1, T0
        if T0 not in T0_L:
            T0_L.append(T0)
        if T0 not in readCntD:
            readCntD[T0] = {}
            baseCntD[T0] = {}
        T0_1 = T0_1.replace(T0, "")
        if not T0_1:
            T0_1 = '1'
        else:
            T0_1 = T0_1[1:]
        #sampleRepD[T0] = sampleRepD.get(T0, 0)+1
        #T0_1 = str(sampleRepD[T0])
        #print >>sys.stderr, T0_1_1, T0_1

        if T0_1 not in T0_1_L:
            T0_1_L.append(T0_1)
        readCntD[T0][T0_1] = sampleL[1]
        baseCntD[T0][T0_1] = str(baseCntD_T0_1[labelMapD[T0_1_1]])
    #--------------------------------------
    #print >>sys.stderr, readCntD
    #print >>sys.stderr, baseCntD

    #sys.exit(1)

    #print tableL
    print "Table: (\#tab:seq-sta-sum) Summary of sequencing reads 测序量总结 \
(对于双端测序,  *\_1* 表示左端reads, *\_2* 表示右端reads) \n"
    print '\n'.join(transferListToMultiLTable(tableL))
    print 

    
    tardir = "1_sequencing_quality_check"
    targetdir = main_dir+'/'+tardir+'/'
    os.system("mkdir -p "+targetdir)

    read_cnt_file = prefix + ".read_cnt.xls"
    read_cnt_fh = open(read_cnt_file, 'w')
    read_cnt_pdf = read_cnt_file + ".barplot.pdf"

    print >>read_cnt_fh, "Sample\t%s" % '\t'.join(T0_1_L)
    for T0 in T0_L:
        print >>read_cnt_fh, "%s\t%s" % (T0, 
            '\t'.join([readCntD[T0].get(i,'0') for i in T0_1_L]))

    read_cnt_fh.close()
    
    if len_labelL <= 60:
        cmd = "s-plot barPlot -f %s -d dodge -y 'Number of sequenced reads' -R 45 -P none -w 30" % read_cnt_file
    else:
        height = len_labelL // 5
        if height < 10:
            height = 10
        cmdL = ["s-plot barPlot -f", read_cnt_file, 
                "-d dodge -y 'Number of sequenced reads' -P none -w 20",
                "-F TRUE -u", str(height)
               ]
        cmd = ' '.join(cmdL)
    os.system(cmd)


    base_cnt_file = prefix + ".base_cnt.xls"
    base_cnt_fh = open(base_cnt_file, 'w')
    base_cnt_pdf = base_cnt_file + ".barplot.pdf"

    print >>base_cnt_fh, "Sample\t%s" % '\t'.join(T0_1_L)
    for T0 in T0_L:
        print >>base_cnt_fh, "%s\t%s" % (T0, 
            '\t'.join([baseCntD[T0].get(i,'0') for i in T0_1_L]))

    base_cnt_fh.close()
    
    if len_labelL <= 60:
        cmd = "s-plot barPlot -f %s -d dodge -y 'Number of sequenced bases' -R 45 -P none -w 30" % base_cnt_file
    else:
        #height = len_labelL // 3
        #if height < 10:
        #    height = 10
        cmdL = ["s-plot barPlot -f", base_cnt_file, 
                "-d dodge -y 'Number of sequenced bases' -P none -w 20",
                "-F TRUE -u", str(height)
               ]
        cmd = ' '.join(cmdL)
    os.system(cmd)

    read_cnt_melt_file = prefix + ".read_cnt_melt.xls"
    read_cnt_melt_fh = open(read_cnt_melt_file, 'w')
    read_cnt_melt_pdf = read_cnt_melt_file + ".densityHist.pdf"

    print >>read_cnt_melt_fh, "variable\tvalue"
    for T0 in T0_L:
        for i in T0_1_L:
            cnt = int(readCntD[T0].get(i, '0'))
            if cnt:
                print >>read_cnt_melt_fh, \
                    "Number of sequenced reads (million)\t%s" % (\
                    str(cnt/(10**6)))

    read_cnt_melt_fh.close()
    cmdL = ["s-plot densityHistPlot -f", read_cnt_melt_file, "-d hist -g ..count..", 
            "-v TRUE -P none -R 45 -x 'Number of sequenced reads (million)'", 
            "-y 'Sample count' -a 1"
            ]
    cmd = ' '.join(cmdL)
    os.system(cmd)

    base_cnt_melt_file = prefix + ".base_cnt_melt.xls"
    base_cnt_melt_fh = open(base_cnt_melt_file, 'w')
    base_cnt_melt_pdf = base_cnt_melt_file + ".densityHist.pdf"

    print >>base_cnt_melt_fh, "variable\tvalue"
    for T0 in T0_L:
        for i in T0_1_L:
            cnt = int(baseCntD[T0].get(i, '0'))
            if cnt:
                print >>base_cnt_melt_fh, \
                    "Number of sequenced bases (G)\t%s" % (\
                    str(cnt/(10**9)))

    base_cnt_melt_fh.close()
    cmdL = ["s-plot densityHistPlot -f", base_cnt_melt_file, "-d hist -g ..count..", 
            "-v TRUE -P none -R 45 -x 'Number of sequenced bases (G)'", 
            "-y 'Sample count' -a 1"
            ]
    cmd = ' '.join(cmdL)
    os.system(cmd)

    copy(targetdir, read_cnt_file, base_cnt_file)
    copypdf(targetdir, read_cnt_pdf, base_cnt_pdf, read_cnt_melt_pdf, base_cnt_melt_pdf)

    #print "## Visualize statistics of sequenced reads\n"

    #knitr_read_txt(main_dir, curation_label)
    
    print "所有样品测序reads数目和碱基数目分布，用于查看是否存在测序量异常的样品 (Figure \@ref(fig:read-num-all-sample) and \@ref(fig:base-num-all-sample))。\n"

    print "\n(ref:read-num-all-sample) Distribution of sequenced reads in all samples. Vertical line represents average of sequenced reads of all samples. 1 Million = 10^6^. \
[PDF](%s/%s)\n" % (tardir, read_cnt_melt_pdf)

    print '''
```{r read-num-all-sample, fig.cap="(ref:read-num-all-sample)"}
knitr::include_graphics("%s/%s")
```

''' % (tardir, read_cnt_melt_pdf.replace('pdf', 'png'))

    print "(ref:base-num-all-sample) Distribution of sequenced bases in all samples. Vertical line represents average of sequenced bases of all samples. 1G = 10^9^. \
[PDF](%s/%s)\n" % (tardir, base_cnt_melt_pdf)

    print '''```{r base-num-all-sample, fig.cap="(ref:base-num-all-sample)"}
knitr::include_graphics("%s/%s")
```

''' % (tardir, base_cnt_melt_pdf.replace('pdf', 'png'))

    print "单个样品测序Reads数目和碱基分布，如果一个样品有多个重复，会以不同颜色的柱子展示 (Figure \@ref(fig:read-num-per-sample))。\n"

    print "(ref:read-num-per-sample) Number of sequenced reads. Color bars represent \
each replicate. 每个样品各个重复测序获得的reads数目。1 Million = 10^6^.\
 [PDF](%s/%s) [SOURCE](%s/%s)\n" % (tardir, read_cnt_pdf, tardir, read_cnt_file)

    print '''```{r read-num-per-sample, fig.cap="(ref:read-num-per-sample)"}
knitr::include_graphics("%s/%s")
```

''' % (tardir, read_cnt_pdf.replace('pdf', 'png'))

    print "(ref:base-num-per-sample) Number of sequenced bases. Color bars represent \
each replicate. 每个样品各个重复测序获得的碱基数。1 G = 10^9^.\
 [PDF](%s/%s) [SOURCE](%s/%s)\n" % (tardir, base_cnt_pdf, tardir, base_cnt_file)

    print '''```{r base-num-per-sample, fig.cap="(ref:base-num-per-sample)"}
knitr::include_graphics("%s/%s")
```

''' % (tardir, base_cnt_pdf.replace('pdf', 'png'))

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


