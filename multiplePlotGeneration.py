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
    This is desigend to generate multiple plots using same command except file name given in JSON file.

Input JSON:

[
    {"grp1": "grp1_matrix"}, 
    {"grp2": "grp2_matrix"}, 
    {"grp3": "grp3_matrix"}
]

Output JSON (as input for ListFileSummaryDoc.py)

{
    "title": "Document title with well number at the beginning to specify the header level", 
    "description": "A paragraph in markdown format to describe this part of analysis.", 
    "tag": "A unique string to label this part of analysis, only alphabets accepted.", 
    "grp1": [
        {"xls": "table file", "pdf": "pdf file", "png":"png_file", "label":"label name", "legend": "A string to describe the pic or table."}, 
        {"xls": "table file", "pdf": "pdf file", "png":"png_file", "label":"label name", "legend": "A string to describe the pic or table."}
        ], 
    "comment": {
        "value": "xls, pdf, png, label, legend are all optional. If only <xls> os given, this file will be directly listed. If <pdf> is given, the png will be generated for shwoing. If <png> is given, the png will be generated for showing. If no <legen>, file name will be used as legend."
    }

}

'''

import sys
import os
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
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="JSON file in format described above.")
    parser.add_option("-t", "--title", dest="title",
        help="Title for this process. Starts with '#'")
    parser.add_option("-T", "--tag", dest="tag",
        help="Tag for this process. Only alphabets allowed.")
    parser.add_option("-d", "--description", dest="description",
        help="Markdown supported description for this process. ")
    parser.add_option("-p", "--plot-command", dest="plot_command",
        help="A s-plot command for plotting with last parameter as filename. Like <s-plot lines -m TRUE -a frag_len -x \"Fragment size (bp)\" -y \"Frequency (%) [lower subplot] or raw counts [upper subplot]\" -f> or simplify <s-plot bars [options] -f>. The needed file name will be directed added to the command to run.")
    parser.add_option("-n", "--no-plot", dest="plot_no",
        default=False, action="store_true", 
        help="Do not do real plot only get JSON file. Usefull when only description words changed.")
    parser.add_option("-l", "--label", dest="plot_label",
        help="As each s-plot command will add a tag to output plot file. Supply the added tag for construct plot filename. Such as <lines>")
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
    file  = options.filein
    title = options.title
    tag   = options.tag
    description  = options.description
    plot_command = options.plot_command
    plot_label   = options.plot_label
    plot_no   = options.plot_no
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    inputL  = json_load(open(file))
    outputD = {}
    outputD['title'] = title
    outputD['tag'] = tag
    outputD['description'] = description
    outputD["value"] = []
    #-------------END reading file----------
    for inputD in inputL:
        for grp, mat in inputD.items():
            plotD = {}
            plotD["xls"] = mat
            plotD["label"] = grp
            plotD["pdf"] = mat+'.'+plot_label+'.pdf'
            cmd = plot_command.strip()+' '+mat
            if not plot_no:
                os.system(cmd)
            outputD['value'].append(plotD)
    ###--------multi-process------------------
    output = file + '.'+plot_label+'.plot.json'
    json_dump(outputD, open(output, 'w'), indent=3)
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


