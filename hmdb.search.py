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
    This is designed to search 'hmdb_metabolites.table' to get the
    KEGG ID of metabolites. Also one can change to get other columns.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from bs4 import BeautifulSoup

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
        metavar="FILEIN", help="One column file with metabolite names.")
    parser.add_option("-f", "--hmdb-table", dest="hmdb",
        metavar="HMDB-TABLE",
        default="/MPATHB/resource/compound/hmdb_metabolites.table",  
        help="Default /MPATHB/resource/compound/hmdb_metabolites.table.")
    parser.add_option("-c", "--column", dest="col",
        default=5, help="1-based number to indicate the column \
one wants to get. Default 5 to get the kegg compound id.")
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
    hmdb = options.hmdb
    col  = int(options.col) - 1
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    metaboliteD = dict([[line.strip().upper(), line.strip()] for line in fh])
    #print metaboliteD
    header = 1
    for line in open(hmdb):
        if header:
            header -= 1
            continue
        #----------------------
        lineL = line.strip().split('\t')
        last_item = lineL.pop(-1)
        lineL.extend(last_item.split('___'))
        for item in lineL:
            item = item.upper()
            if item in metaboliteD:
                print "%s\t%s" % (metaboliteD[item], lineL[col])
                metaboliteD.pop(item)
        #-----------------------------------
        if not metaboliteD:
            break
    #-------------END reading file----------
    if metaboliteD:
        for key, item in metaboliteD.items():
            print "%s\t%s" % (item, 'Not found')

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


