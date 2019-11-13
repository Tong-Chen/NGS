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

JSON format:

{
    "title": "Document title with well number at the beginning to specify the header level", 
    "description": "A paragraph in markdown format to describe this part of analysis.", 
    "tag": "A unique string to label this part of analysis, only alphabets accepted.", 
    "value": [
        {"xls": "table file", "pdf": "pdf file", "png":"png_file", "label":"sub title", "other_key":"other_type_file", "legend": "A string to describe the pic or table."}, 
        {"xls": "table file", "pdf": "pdf file", "png":"png_file", "label":"sub title", "other_key":"other_type_file", "legend": "A string to describe the pic or table."}
        ], 
    "comment": {
        "value": "xls, pdf, png, legend are all optional. If only <xls> or <other_key (other_key can be any string)> is given, this file will be directly listed. If <pdf> is given, the png will be generated for showing. The <xls> will be treated as source file for <pdf>. If <png> is given, the png will be generated for showing. The <xls> will be treated as source file for <png>. If no <legend>, file name will be used as legend."
    }

}
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
    parser.add_option("-f", "--json-file", dest="filein",
        metavar="FILEIN", help="JSON file containing data information.")
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

def generateDoc(report_dir, report_sub_dir, fileD):
    dest_dir = report_dir+'/'+report_sub_dir+'/'
    os.system('mkdir -p '+dest_dir)

    print "{}\n".format(fileD['title'])
    level = fileD['title'].find(' ')
    if level == -1:
        level = title.count('#')
    if level < 2:
        level = 2
    #curation_label = "Sequencing_saturation_estimation"
    #knitr_read_txt(report_dir,  curation_label)
    
    print """
{}

""".format(fileD["description"])

    tag = fileD["tag"]
    
    valueL = fileD["value"]
    
    count = 0
    for aDict in valueL:
        count += 1
        specialL = ['pdf', 'png', 'xls', 'label']
        other_keyL = [i for i in aDict.keys() if i not in specialL]
        label = aDict.get('label', '')
        if label:
            print 
            print '#'*level, label
            print
        pdf = aDict.get('pdf', '')
        png = aDict.get('png', '')
        if pdf:
            copypdf(dest_dir, pdf)
            pdf_name = os.path.split(pdf)[-1]
            pdf = report_sub_dir+'/'+pdf_name
            png = pdf.replace('pdf', 'png')
        else:
            pdf = '#'
        if png:
            png_name = os.path.split(png)[-1]
            if pdf == '#':
                copypng(dest_dir, png)
                png = report_sub_dir+'/'+png_name
        xls = aDict.get('xls', '')
        if xls:
            copy(dest_dir, xls)
            xls_name = os.path.split(xls)[-1]
            xls = report_sub_dir+'/'+xls_name
        else:
            xls = '#'
        legend = aDict.get('legend', '')
        
        if png:
            if not legend:
                legend = png_name.replace('.png', '') 
            print "(ref:{tag}-{count}) {legend} ([XLS table]({xls}); [PDF pic]({pdf}))\n".format(tag=tag, count=count, legend=legend, xls=xls, pdf=pdf)

            print '''```{{r {tag}-{count}, fig.cap="(ref:{tag}-{count})"}}
knitr::include_graphics("{png}")
```
'''.format(tag=tag, count=count, png=png)
        elif xls not in ['', '#']:
            if not legend:
                legend = xls_name
            print "* [{legend}]({xls})\n".format(legend=legend, xls=xls)
        if other_keyL:
            for other_key in other_keyL:
                legend = aDict.get('legend', '')
                other_file = aDict[other_key]
                copy(dest_dir, other_file)
                other_file_name = os.path.split(other_file)[-1]
                other_file = report_sub_dir+'/'+other_file_name
                if not legend:
                    legend = other_file_name
                print "* [{legend}]({other_file})\n".format(legend=legend, other_file=other_file)
        #------------------------------------------
#--------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    verbose = options.verbose
    report_dir = options.report_dir
    report_sub_dir = options.report_sub_dir
    global debug
    debug = options.debug
    
    fileD = json_load(open(file))
    generateDoc(report_dir, report_sub_dir, fileD)
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


