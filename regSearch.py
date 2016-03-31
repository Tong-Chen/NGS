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
    This is designed to do regular expression pattern search in FASTA
    files.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from bs4 import BeautifulSoup
import re

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
    usages = "%prog -i file -p 'CW.{2,5}D.DD[A-Z][0-9]\\s'"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A FASTA file.")
    parser.add_option("-p", "--pattern", dest="pat",
        metavar="PATTERN", help="The regular expression pattern.")
    parser.add_option("-m", "--match", dest="match",
        default=1, help="Output matched sequences (Default). \
Accept 0 to output unmatched sequences.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def patSearch(pat, sequence):
    pat_match = pat.findall(sequence)
    return ';'.join(pat_match)
#---END patSearch---------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    pat  = options.pat
    output_match = int(options.match)
    pat  = re.compile(r"%s" % pat)
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    key = ''
    for line in fh:
        if line[0] == '>':
            if key:
                seq = ''.join(seqL)
                matched = patSearch(pat, seq)
                if output_match and matched:
                    print "%s|%s\n%s" % (key, matched, seq)
                elif (not output_match) and (not matched):
                    print "%s\n%s" % (key, seq)
                #--------------------------------------------
            key = line.strip()
            seqL = []
        else:
            seqL.append(line.strip())
    #-------------END reading file----------
    if key:
        seq = ''.join(seqL)
        matched = patSearch(pat, seq)
        if output_match and matched:
            print "%s|%s\n%s" % (key, matched, seq)
        elif (not output_match) and (not matched):
            print "%s\n%s" % (key, seq)
        #--------------------------------------------
    #--------------------------------------------
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


