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
    This is designed to count the number of CHH, CHG, CpG (H
    represents A, T, C) of one genome or other FASTA sequences.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import re
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
        metavar="FILEIN", help="FASTA sequence file")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def count_c(sequence, chg_p, chh_p, chg_p2, chh_p2):
    c   = sequence.count('C')
    g   = sequence.count('G')
    cg  = sequence.count('CG')
    chg = len(chg_p.findall(sequence))
    chh = len(chh_p.findall(sequence))
    chg_2 = len(chg_p2.findall(sequence))
    chh_2 = len(chh_p2.findall(sequence))
    return c, g, cg, chg, chh, chg_2, chh_2
#-------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    chg_p = re.compile(r'C[ACT]G')
    chh_p = re.compile(r'C[ACT][ACT]')
    chg_p2 = re.compile(r'C[AGT]G')
    chh_p2 = re.compile(r'[AGT][AGT]G')
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    countL = []
    seqL = []
    for line in fh:
        if line[0] == '>' and seqL:
            seq = ''.join(seqL)
            c, g, cg, chg, chh, chg2, chh2 \
                = count_c(seq, chg_p, chh_p, chg_p2, chh_p2)
            countL.append([c, g, cg, chg, chh, chg2, chh2])
            seqL = []
        #---------------------------------
        else:
            seqL.append(line.strip())
    #-------------END reading file----------
    if seqL:
        seq = ''.join(seqL)
        c, g, cg, chg, chh, chg2, chh2 \
            = count_c(seq, chg_p, chh_p, chg_p2, chh_p2)
        countL.append([c, g, cg, chg, chh, chg2, chh2])
    
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    C = G = CG = CHG = CHH = CHG2 = CHH2 = 0
    for c, g, cg, chg, chh, chg2, chh2 in countL:
        C  += c
        G  += g
        CG += cg
        CHG += chg
        CHH += chh
        CHG2 += chg2
        CHH2 += chh2
    #----------------------------------------
    print "C:%d" % C
    print "G:%d" % G
    print "CpG:%d" % CG
    print "CHG(+):%d" % CHG
    print "CHG(-):%d" % CHG2
    print "CHH(+):%d" % CHH
    print "CHH(-):%d" % CHH2

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


