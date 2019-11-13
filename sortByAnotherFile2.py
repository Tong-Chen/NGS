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
    This is designed to sort one file by the other file.



classCode

.
u
r
p
i
s
x
o
c
e
j
=

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

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
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="File waiting for sort")
    parser.add_option("-c", "--sort-col", dest="sort_col",
        type="int", default=1, help="Specify the column to sort on. Default 1 meaning sort file on the first column.")
    parser.add_option("-n", "--number-head", dest="num_head",
        type="int", default=1, help="Specify the number of headlines to skip.")
    parser.add_option("-o", "--order-file", dest="order_file",
        help="The file containing ordered lists (one in each line).")
    parser.add_option("-C", "--order-col", dest="order_col",
        type="int", default=1, help="Specify the column in <order_file> used for sort <input-file>. Default 1 meaning used the first column.")
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
    file       = options.filein
    sort_col   = options.sort_col-1
    header     = options.num_head
    order_file = options.order_file
    order_col  = options.order_col-1
    verbose    = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    orderL = [line.strip().split()[order_col] for line in open(order_file)]
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    fileD = {}
    for line in fh:
        line = line.strip()
        if header:
            print line
            header -= 1
            continue
        key = line.split()[sort_col]
        if key not in fileD:
            fileD[key] = []
        fileD[key].append(line)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    for key in orderL:
        try:
            tmpL = fileD.pop(key)
        except KeyError:
            continue
        print '\n'.join(tmpL)
    assert len(fileD) == 0, fileD.keys()
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


