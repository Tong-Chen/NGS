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
    This is used to integrate sample info into FASTA seq names.

Suppose there is a file <SRR037890.fasta>, which contents are:

>comp0_c0_seq1 len=597 path=[4541:0-103 5688:104-596]
ACGCAGAGTACGCGGGGGAGATGGGATGGGACTCAAAGCGCCTGTGGTCTACATTATGCT
TCACTCTCTATACTTATAAAAACAGTACAATCTTCAAATCCTTTGTTCTCATTCCATATC
ATGGGCTCACACTATCCATGAGCCCGAGAAAGGGTGTTTACGTCATTGCCTACAACG
>comp1_c0_seq1 len=509 path=[5:0-269 1771:270-293 9016:294-508] long_read_mappings: {PairPath [_paths=[[5, 1771, 9016], []]]=[SRR037890.141645/F]}
AACGCAGAGTACGCGGGGACCAGACACTTTCTCCTTCTCTGTTTATCAAACCCACCTTGA
TTATATTAGTAAACAAAAAAAAAAAAAAA
>comp1_c0_seq2 len=373 path=[6337:0-348 1771:349-372]
TGAAGATGCATTTGATCCAAATCGTGTGAAATGTGTGGTGGATAATCGTGGCTATGCAAT
GGATTCAGAGCTATGATTCAAAGTTCCCAGAATATATTCAGAGCTGCCACCTACCCC

The output would be 

>SRR037890@comp0_c0_seq1
ACGCAGAGTACGCGGGGGAGATGGGATGGGACTCAAAGCGCCTGTGGTCTACATTATGCT
TCACTCTCTATACTTATAAAAACAGTACAATCTTCAAATCCTTTGTTCTCATTCCATATC
ATGGGCTCACACTATCCATGAGCCCGAGAAAGGGTGTTTACGTCATTGCCTACAACG
>SRR037890@comp1_c0_seq1 len=509 path=[5:0-269 1771:270-293 9016:294-508] long_read_mappings: {PairPath [_paths=[[5, 1771, 9016], []]]=[SRR037890.141645/F]}
AACGCAGAGTACGCGGGGACCAGACACTTTCTCCTTCTCTGTTTATCAAACCCACCTTGA
TTATATTAGTAAACAAAAAAAAAAAAAAA
>SRR037890@comp1_c0_seq2 len=373 path=[6337:0-348 1771:349-372]
TGAAGATGCATTTGATCCAAATCGTGTGAAATGTGTGGTGGATAATCGTGGCTATGCAAT
GGATTCAGAGCTATGATTCAAAGTTCCCAGAATATATTCAGAGCTGCCACCTACCCC

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
        metavar="FILEIN", help="A FASTA file.")
    parser.add_option("-l", "--label", dest="label",
        help="Optional. Default file name")
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
    label = options.label
    if not label:
        label = file.rsplit('.',1)[0]
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
        if line[0] == '>':
            id = line.split()[0][1:]
            print ">{}___{}".format(label, id)
        else:
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

