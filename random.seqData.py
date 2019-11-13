#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2017, 陈同 (chentong_biology@163.com).  
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
    parser.add_option("-s", "--samplecnt", dest="filein",
        type="int", metavar="FILEIN", help="Sample conut")
    parser.add_option("-n", "--number", dest="number",
        default=20100100,
        type="int", help="Expected sequence number")
    parser.add_option("-l", "--length", dest="length",
        default=150,
        type="int", help="Expected read length")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

from random import randint

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    sampleCnt = options.filein
    expectedSeq = options.number
    minSeqnum = int(expectedSeq * (1-0.1**8)+0.5)
    maxSeqNum = int(expectedSeq * (1+0.1**3)+0.5)
    expectedLen = options.length
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    print "Sample\tRaw reads\tClean reads\tRaw bases\tClean bases\tQ20 (%)\tQ30 (%)\tGC content (%)"
    for i in range(sampleCnt):
        sampleid = "%s%03d" % ("sample", i+1)
        raw_reads = randint(minSeqnum, maxSeqNum)
        rand_ratio = randint(9999980,10000000) / 10000000.0
        clean_reads = int(raw_reads * rand_ratio)
        raw_base = "%.4f" % (raw_reads * expectedLen * 2.0 / (10 ** 9))
        clean_base = "%.4f" % (clean_reads * expectedLen * 2.0 * rand_ratio / (10 ** 9))
        Q20 = "%.2f" % (randint(9600,9710)/100.0)
        Q30 = "%.2f" % (randint(9200,9400)/100.0)
        GC = "%.2f" % (randint(4400,4800)/100.0)
        print '\t'.join([sampleid, str(raw_reads),str(clean_reads), str(raw_base), str(clean_base), 
            str(Q20),str(Q30), str(GC)])
    #-----------end close fh-----------
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
    if verbose:
        print >>sys.stderr,            "--Successful %s" % strftime(timeformat, localtime())

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " %         (' '.join(sys.argv), startTime, endTime)
    fh.close()
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------
#

