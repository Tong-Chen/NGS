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
    This is designed to transfer the output of eprimer32 to input of primersearch.

Input:

# EPRIMER32 RESULTS FOR comp31_c0_seq1

#                      Start  Len   Tm     GC%   Sequence

   1 PRODUCT SIZE: 203
     FORWARD PRIMER     151   20  50.29  30.00  ATCATGCTTTTCCTATGAAA

     REVERSE PRIMER     334   20  50.29  30.00  AAAAACAATTCCCATTCAAG


   2 PRODUCT SIZE: 204
     FORWARD PRIMER     150   20  50.29  30.00  AATCATGCTTTTCCTATGAA

     REVERSE PRIMER     334   20  50.29  30.00  AAAAACAATTCCCATTCAAG


Output


# This is my primer file
D1S243  cacacaggctcacatgcc      gctccagcgtcatggact
D1S468 aattaaccgttttggtcct     gcgacacacacttccc 
D1S2845 ccaaagggtgcttctc        gtggcattccaacctc
D1S1608 gatggcttttggggactatt    cactgagccaagtgacacag
D1S2893 aaaacatcaactctcccctg    ctcaaaccccaataagcctt
D1S2660 cacacatgcacatgcac       agtgacaccagcaggg



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
        metavar="FILEIN", help="EPRIMER32 output file")
    parser.add_option("-o", "--outputfile", dest="outputfile",
        help="Output file")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def transferPrimer3OutputToPrimersearchInput(primer3, primersearch):
    fh = open(primer3)
    fh_out = open(primersearch,'w')
    line = fh.readline()
    line = fh.readline()
    assert line.find("EPRIMER32 RESULTS FOR") != -1, "Wrong format" + primer3
    seq_name = line.strip()[24:]
    
    #primerL = []
    count = 1
    for line in fh:
        if line.find("FORWARD PRIMER") != -1:
            forward = line.strip().split()[-1]
        if line.find("REVERSE PRIMER") != -1:
            reverse = line.strip().split()[-1]
            #tmpL = [seq_name+'@'+str(count), forward, reverse]
            print >>fh_out, "{}@{}\t{}\t{}".format(seq_name,count,forward, reverse)
            #primerL.add(tmpL)
            count += 1
    #------------------------------------------
    fh.close()
    fh_out.close()
#---------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    output = options.outputfile
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    transferPrimer3OutputToPrimersearchInput(file,output)
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

