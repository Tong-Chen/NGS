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

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
import gzip

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
        metavar="FILEIN", help="FASTQ prefix")
    parser.add_option("-t", "--type", dest="type",
        help="PE or SE")
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
    file = options.filein
    type = options.type
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if type == 'SE':
        file = file + '.fq.gz'
        i = 0
        for line in gzip.open(file, 'rb'):
            i += 1
            if i % 4 == 1:
                key = line
                info = key +'...'+ str(i)
                assert line[0] == '@', info 
            elif i % 4 == 2:
                assert line[0] in ['A', 'C', 'G', 'T', 'N'], info
            elif i % 4 == 3:
                assert line[0] == '+', info
    if type == 'PE':
        file1 = file +'_1.fq.gz'
        file2 = file +'_2.fq.gz'
        
        file1_fh = gzip.open(file1, 'rb')
        file2_fh = gzip.open(file2, 'rb')
        i = 0
        for line1 in file1_fh:
            i += 1
            line2 = file2_fh.readline()
            if i % 4 == 1:
                key1 = line1.split()[0]
                key2 = line2.split()[0]
                assert key1 == key2, key1+'...'+key2+'...'+str(i)
                info = key1 + '...' + str(i)
                assert line1[0] == '@', info 
                assert line2[0] == '@', info 
            elif i % 4 == 2:
                assert line1[0] in ['A', 'C', 'G', 'T', 'N'], info
                assert line2[0] in ['A', 'C', 'G', 'T', 'N'], info
            elif i % 4 == 3:
                assert line1[0] == '+', info
                assert line2[0] == '+', info

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


