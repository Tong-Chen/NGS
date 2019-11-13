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
    This is designed to extract GTF attributes including mRNA length, CDS length, exon counts,  exons size, intron size from BED12 file.
'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
import re

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
        metavar="FILEIN", help="<,> or < > separated file lists")
    parser.add_option("-l", "--labels", dest="labels",
        metavar="FILEIN", help="<,> or < > separated file labels")
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
    fileL = re.split(r"[, ]*", options.filein)
    labelL = re.split(r"[, ]*", options.labels)
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    print "type\tvalue\tvariable"
    for file, label in zip(fileL, labelL):
        for line in open(file):
            lineL = line.strip().split('\t')
            gene_len = int(lineL[2])-int(lineL[1])
            #CDS_len  = int(lineL[7])-int(lineL[6])
            exon_count = int(lineL[9])
            exonSizeL = [int(i) for i in lineL[10].strip(',').split(',')]
            mRNA_len = sum(exonSizeL)
            exonStartL = [int(i) for i in lineL[11].strip(',').split(',')]
            intronSizeL = []
            for i in range(1, exon_count):
                intronSizeL.append(exonStartL[i]-exonSizeL[i-1]-exonStartL[i-1])
            print "Gene_len\t{}\t{}".format(gene_len, label)
            print "mRNA_len\t{}\t{}".format(mRNA_len, label)
            print "exon_count\t{}\t{}".format(exon_count, label)
            for exon_size in exonSizeL:
                print "exon_size\t{}\t{}".format(exon_size, label)
            for intron_size in intronSizeL:
                print "intron_size\t{}\t{}".format(intron_size, label)
            
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


