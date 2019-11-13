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
    This is designed to generate a json file for given lists of files and the json file will be used by <ListFileSummaryDoc.py>.

The output json would be in format like:


JSON format:

{
    "title": "Document title with well number at the beginning to specify the header level", 
    "description": "A paragraph in markdown format to describe this part of analysis.", 
    "tag": "A unique string to label this part of analysis, only alphabets accepted.", 
    "value": [
        {"xls": "table file", "pdf": "pdf file", "png":"png_file", "other_key":"other_type_file", "legend": "A string to describe the pic or table."}, 
        {"xls": "table file", "pdf": "pdf file", "png":"png_file", "other_key":"other_type_file", "legend": "A string to describe the pic or table."}
        ], 
    "comment": {
        "value": "xls, pdf, png, legend are all optional. If only <xls> or <other_key (other_key can be any string)> is given, this file will be directly listed. If <pdf> is given, the png will be generated for showing. The <xls> will be treated as source file for <pdf>. If <png> is given, the png will be generated for showing. The <xls> will be treated as source file for <png>. If no <legend>, file name will be used as legend."
    }

}

'''

import sys
import os
import json
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0


def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-string", dest="filein",
        metavar="FILEIN", 
        help="Input file can be in format like <xls:a.xls,pdf:a.pdf,png:a.png,bw:a.bw,legend:a.legend (no , or : or ; allowed in legend);xls:b.xls,pdf:b.pdf;bw:c.bw,legend:c.legend>. <;> separates each sample which will be saved in one independent dict. <,> separates each file. <:> separate file type and file name which will be treated as a key-value pair.")
    #parser.add_option("-t", "--type", dest="type",
    #    metavar="FILEIN", help="Type for each part of input file like <xls> or <xls;pdf> relative to different input format.")
    parser.add_option("-T", "--title", dest="title",
        help="Document title with well number at the beginning to specify the header level. Currently only first level header supported.")
    parser.add_option("-d", "--description", 
        dest="description", 
        help="A paragraph in markdown format to describe this part of analysis.")
    parser.add_option("-t", "--tag", dest="tag", 
        help="A unique string to label this part of analysis,  only alphabets accepted.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    string = options.filein
    title  = options.title
    description = options.description
    tag = options.tag
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    aDict = {}
    aDict['title'] = title
    aDict['description'] = description
    aDict['tag'] = tag
    aDict['value'] = []

    eachItemL = string.split(';')
    for eachItem in eachItemL:
        if not eachItemL:
            continue
        tmpL = [i.strip().split(':') for i in eachItem.split(',') if i.strip()]
        tmpD = dict(tmpL)
        aDict['value'].append(tmpD)
    #-----------------------------------
    print json.dumps(aDict, indent=3)
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


