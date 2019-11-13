#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright 2018, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    this is designed to extract go id (alt_id), go term, go category information from go.obo file download from http://purl.obolibrary.org/obo/go.obo. 
'''

import sys
import os
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
        metavar="FILEIN", help="go.obo")
    #parser.add_option("-n", "--number", dest="number",
    #    type="int", help="Supply an int number")
    #parser.add_option("-c", "--choice", dest="choice",
    #    type="choice", choices=["a", "b", "c"], 
    #    default="a", help="Supply an int number")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
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
        if line.startswith('id: '):
            id1 = line.strip().split(": ")[1]
            name = ''
            namespace = ''
        elif line.startswith('name: '):
            name = line.strip().split(": ")[1]
        elif line.startswith('namespace: '):
            namespace = line.strip().split(": ")[1]
            assert id1 and name and namespace, id1 + name + namespace
            print "\t".join([id1, name, namespace])
        elif line.startswith('alt_id: '):
            alt_id = line.strip().split(": ")[1]
            assert alt_id and name and namespace, alt_id + name + namespace
            print "\t".join([alt_id, name, namespace])

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
        #print("--Successful %s" % strftime(timeformat, localtime()), file=sys.stderr)
        print >>sys.stderr, "--Successful %s" % strftime(timeformat, localtime())

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    with open('python.log', 'a') as fh:
        #print ("%s\n\tRun time : %s - %s " % \
        #(' '.join(sys.argv), startTime, endTime), file=fh)
        print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


