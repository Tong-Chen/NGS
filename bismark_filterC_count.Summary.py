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
    This is designed to summarize results output of files in following format.

    variable    value   sample
    CpG 100 A
    CHG 100 A
    CHH 100 A


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
        metavar="FILEIN", help="`,` or ` ` separated a list of files. ")
    parser.add_option("-l", "--labels", dest="label",
        metavar="LABEL", help="`,` or ` ` separated a list of labels to label each file. It must have same order as files.")
    parser.add_option("-c", "--reads-support", dest="reads_support",
        type='int', default=5, help="Specify the a number of reads to indicate trustable coverage.")
    parser.add_option("-o", "--output-file", dest="out_file",
        help="The name of output files. ")
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


#---------------------------------------
def plot(op):
    cmd = ["s-plot barPlot -f", op, "-m TRUE -a Sample -R 45", 
           "-x 'Samples' -y 'Number of covered cytosine sites'", 
           "-P none -B variable -O 1 -b \"'CpG','CHG','CHH'\"", 
           "-k free_y -u 20"]
    #print >>sys.stderr, ' '.join(cmd)
    os.system(' '.join(cmd))
#-----------------------------------------------------

def mergeFile(fileL, op):
    fh = open(op, 'w')
    print >>fh, "variable\tvalue\tSample"
    for file in fileL:
        header = 1
        for line in open(file):
            if header:
                header -= 1
                continue
            #------------------------
            print >>fh, line, 
    fh.close()
#-----------------------------------

def generateDoc(report_dir, report_sub_dir, op, fileL, labelL, 
        curation_label, reads_support):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)

    fileL, labelL = skipUnExistedOrEmptyFile(fileL[:], labelL[:])
    mergeFile(fileL, op)
    plot(op)
    pdf = op + '.stackBars.pdf'
    copy(dest_dir, op)
    copypdf(dest_dir, pdf)

    print "\n## 甲基化位点数目统计\n"

    knitr_read_txt(report_dir,  curation_label)
    
    print """为了保证结果的可信度，选取有不少于{} reads覆盖的位点用于后续分析，每个样品中检测到的甲基化位点数目如图所示。

""".format(reads_support)

    print "(ref:mc-count-{read}) Summary of detected cytosine sites in each sample. Only cytosine sites with no less than {read} reads support were counted.\n".format(read=reads_support)
        
    png = report_sub_dir+'/'+os.path.split(pdf)[-1].replace('pdf', 'png')

    print '''```{{r mc-count-{read}, fig.cap="(ref:mc-count-{read})"}}
knitr::include_graphics("{png}")
```

'''.format(read=reads_support, png=png)

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
    op = options.out_file
    reads_support = options.reads_support
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
    curation_label = os.path.split(sys.argv[0])[-1].replace('.', '_')
    if doc_only:
        generateDoc(report_dir, report_sub_dir, op, fileL[:], labelL[:], curation_label, 
                reads_support)
        return 0

    generateDoc(report_dir, report_sub_dir, op, fileL[:], labelL[:], curation_label, 
            reads_support)
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


