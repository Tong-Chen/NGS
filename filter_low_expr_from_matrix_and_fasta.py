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
    This is designed to filter out low expressed genes from the
    expression matrix.
    For Trinity assembled transcripts, the FASTA sequence can also be
    filtered according to matrix.

expr_matrix for -i

Gene    Samp1   Samp2   Samp3
a   0   1   2
b   0   1   2
c   0   1   2
d   0   1   2

Fasta file given to -f (optional)
# The original Fasta file will be saved and filtered fasta will be
# output to original file.

# Only strings between the beginning > and the first blank 
# will be treated as sequence names.

# The output will keep the original format with no information
# discarded for each sequence.

>a  len=563 path=[541:0-562] [-1,  541,  -2]
ACGATCGAGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
ACGATCGAGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>B  len=563 path=[541:0-562] [-1,  541,  -2]
ACGATCGAGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
ACGATCGAGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>d  len=563 path=[541:0-562] [-1,  541,  -2]
ACGATCGAGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
ACGATCAGTCAGGAGTAGGATGTAGGAGTAGCATG

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
        metavar="FILEIN", help="An expression matrix. \
Filtered matrix will be output to STDOUT.")
    parser.add_option("-f", "--fasta-file", dest="fasta",
        help="A FASTA file with formats specified above. \
Filtered FASTA will be saved in original file. \
Original FASTA will be renamed.")
    parser.add_option("-m", "--minimum-expr-value", dest="min_expr",
        default=0, help="The minimum expression level for one \
transcript to be kept.")
    parser.add_option("-l", "--no-less-or-larger", dest="no_less",
        default=1, help="Transcripts with expression in at least one \
sample no less than the minimum expression level will be kept by \
default. Given 0 here to exclude transcrpts with expression value \
exactly the same as given <minmum-expr-value>.")
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
    fasta = options.fasta
    min_expr = float(options.min_expr)
    no_less = int(options.no_less)
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    header = 1
    os.system("/bin/mv -f %s %s.bak" % (file, file))
    matrix_fh = open(file, 'w')
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file+'.bak')
    #--------------------------------
    savedD = {}
    for line in fh:
        if header:
            print >>matrix_fh, line,
            header -= 1
            continue
        #-------------------------------------
        lineL = line.strip().split('\t')
        if no_less:
            tmpL = [1 for i in lineL[2:] if float(i)>=min_expr]
        else:
            tmpL = [1 for i in lineL[2:] if float(i)>min_expr]
        #-----------------------------------
        if tmpL:
            print >>matrix_fh, line,
            savedD[lineL[0]] = 1
        #----------------------
    #-------------END reading file----------
    #print savedD
    if fasta:
        os.system("/bin/mv -f %s %s.bak" % (fasta, fasta))
        fasta_fh = open(fasta, 'w')
        for line in open(fasta+'.bak'):
            if line[0] == '>':
                keep = 0
                name = line.split()[0][1:]
                #print name
                if name in savedD:
                    keep = 1
                    print >>fasta_fh, line,
            else:
                if keep:
                    print >>fasta_fh, line,
        #-----------------------------
        fasta_fh.close()
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


