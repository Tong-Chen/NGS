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
    This is designed to summarize results of JSON format.

'''

import sys
import os
from json import dumps as json_dumps
from json import load as json_load
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
    parser.add_option("-p", "--prefix", dest="prefix",
        metavar="prefix", help="Used to set the prefix and get the <summary/prefix.Methylrate_distrib.bp.xls> series file.")
    parser.add_option("-c", "--CpG", dest="CpGIsland",
        default="CpG_islands,CpG_shelves,CpG_shores", 
        help="Used to get <sumamry/prefix.CpG_shelves.CpGisland_Methylrate_distrib_sum.bp.xls> series files. Default <CpG_islands, CpG_shelves, CpG_shores>.")
    parser.add_option("-l", "--rmsk", dest="rmsk",
        default="LINE,SINE,LTR,Retroposon", 
        help="Used to get <sumamry/prefix.LTR.Repeat_Methylrate_distrib_sum.bp.xls> series files. Default <LINE, SINE, LTR, Retroposon>.")
    parser.add_option("-i", "--ICR", dest="ICR",
        default="ICRmethylRate", 
        help="Used to get <sumamry/prefix.ICRmethylRate_sum.xls> series files. Default <ICRmethylRate>.")
    parser.add_option("-P", "--promoter", dest="promoter",
        default="promoterMethylRate", 
        help="[Uppercase P] Used to get <sumamry/prefix.promoterMethylRate_sum.xls> series files. Default <ICRmethylRate>.")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='6_customized_result', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '6_customized_result'.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def generateDoc(report_dir, report_sub_dir, prefix, CpGIsland, rmsk, ICR, promoter):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)

    print "\n## 甲基化率分布统计\n"
    
    summary_dir = "summary"

    curation_label = "Methylrate_distribution"
    knitr_read_txt(report_dir,  curation_label)
    
    print "\n### 整体甲基化率分布\n"
    tag = 'methylRateTotal'
    count = 0

    xls = summary_dir+'/'+prefix+'.Methylrate_distrib.bp.xls'
    bp_pdf = xls + '.boxplot.value_Individual.pdf'
    copypdf(dest_dir, pdf)
    bp_pdf = getRelativeDir(bp_pdf, report_sub_dir)
    bp_png = bp_pdf.replace('pdf', 'png')
    
    legend = "CpG methylation rate distribution of all samples."

    print "(ref:{tag}-{count}) {legend} ([XLS table]({xls}); [PDF pic]({pdf}))\n".format(tag=tag, count=count, legend=legend, xls=xls, pdf=bp_pdf)

    print '''```{{r {tag}-{count}, fig.cap="(ref:{tag}-{count})"}}
knitr::include_graphics("{png}")
```
'''.format(tag=tag, count=count, png=bp_png)
    
    if CpGIsland:
        print """ 
# CpG island methylation rate distribution



"""
        CpGIslandL = re.split(r'[, ]', CpGIsland)

    if rmsk
    if ICR
    if promoter):

        #------------------------------------------
#--------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    prefix = options.prefix
    CpGIsland = options.CpGIsland
    rmsk = options.rmsk
    ICR = options.ICR
    promoter = options.promoter

    verbose = options.verbose
    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    global debug
    debug = options.debug
    
    fileD = json_load(open(file))
    generateDoc(report_dir, report_sub_dir, prefix, CpGIsland, rmsk, ICR, promoter)
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


