#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2013, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Functional description:
    This is designed to extract rows containing given words.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from bs4 import BeautifulSoup
import re

#reload(sys)
#sys.setdefaultencoding('utf8')

#from multiprocessing.dummy import Pool as ThreadPool

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
        metavar="FILEIN", help="A matrix file normally with \
one header line.")
    parser.add_option("-p", "--pattern", dest="pattern",
        help="Words or regular expressions used as filter.")
    parser.add_option("-H", "--header", dest="header",
        default=1, help="Output the header line. \
Default the first line will be treated as header line and \
will be output. Accept 0 to ignor header lines.")
    parser.add_option("-c", "--case-insensitive", dest="nocase",
        default=1, help="Default ignore case. \
Accept 0 to be case sensitive.")
    parser.add_option("-f", "--filter-list", dest="filter_l",
        help="Only IDs exist in this file will be output \
once given. Only the first column will be used.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    pat = options.pattern
    header = int(options.header)
    nocase = int(options.nocase)
    filter_l = options.filter_l
    if filter_l:
        filterD = dict([(line.split()[0], 1) for line in open(filter_l)])
    else:
        filterD = ''
    #--------------------------
    if nocase:
        pat_re = re.compile(r'%s' % pat, flags=re.IGNORECASE)
    else:
        pat_re = re.compile(r'%s' % pat)
    #-------------------------------
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    for line in fh:
        if header:
            print line,
            header -= 1
            continue
        #---------------------------------
        if filterD:
            key = line.split()[0]
            if key not in filterD:
                continue
        #--------------------------------
        if pat_re.search(line):
            print line,
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
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


