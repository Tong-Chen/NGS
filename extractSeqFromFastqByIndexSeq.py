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

    index.table

    GATCAG  113F
    TAGCTT  113M

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import gzip
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
        metavar="FILEIN", help="Input fastq file name")
    parser.add_option("-t", "--seq-type", dest="seq_type",
        help="<PE> or <SE>, for <PE> <_1.fq.gz> <_2.fq.gz> would be added to generate full fastq file name. For <SE>, <.fq.gz> would be added to make full fastq file name.")
    parser.add_option("-I", "--index-file", dest="index_file",
        help="A matrix with first column as index sequences and second column as output file name prefix (see above for an example).")
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
    seq_type = options.seq_type
    index_file = options.index_file
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    mapD = {'A':'T', 'G':'C', 'T':'A','C':'G','a':'t','g':'c','t':'a','c':'g'}
    #indexL = [line.split() for line in open(index_file)]
    indexD = {}
    for line in open(index_file):
        seq, name = line.split()
        revseq = [mapD[i] for i in seq]
        revseq.reverse()
        revseq = ''.join(revseq)
        assert seq not in indexD, "Duplicate "+seq
        indexD[seq] = {}
        assert revseq not in indexD, "Duplicate "+revseq
        indexD[revseq] = {}
        if seq_type == "PE":
            tmp_file_1 = gzip.open(name + '_1.fq.gz', 'wb')
            tmp_file_2 = gzip.open(name + '_2.fq.gz', 'wb')
            indexD[seq]['1'] = tmp_file_1
            indexD[seq]['2'] = tmp_file_2
            indexD[revseq]['1'] = tmp_file_1
            indexD[revseq]['2'] = tmp_file_2
        else:
            tmp_file = gzip.open(name + '.fq.gz', 'wb')
            indexD[seq]['1'] = tmp_file
            indexD[revseq]['1'] = tmp_file
    if debug:
        print >>sys.stderr, indexD
    all_index_seqL = indexD.keys()
    if seq_type == "PE":
        inputD = {'1': file+'_1.fq.gz', '2': file+'_2.fq.gz'}
    else:
        inputD = {'1': file+'.fq.gz'}
    #------------------------------------------------------
    if debug:
        print >>sys.stderr, inputD
    for end, file in inputD.items():
        output = 0
        for line in gzip.open(file, 'rb'):
            if line[0] == '@':
                output = 0
                for tmp_seq in all_index_seqL:
                    if line.find(tmp_seq) != -1:
                        output = 1
                        print >>indexD[tmp_seq][end], line,
                        break
            elif output:
                print >>indexD[tmp_seq][end], line,
        #--------------------------------------------
    #0000000000000000000-----------------------
    closedFh = []
    for innerD in indexD.values():
        for fh in innerD.values():
            if fh not in closedFh:
                closedFh.append(fh)
                if debug:
                    print >>sys.stderr, fh
                fh.close()
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


