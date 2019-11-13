#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright 2018, 陈同 (chentong_biology@163.com).  
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
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
from subprocess32 import check_output, PIPE
#from multiprocessing.dummy import Pool as ThreadPool

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
        metavar="FILEIN", help="FASTA file")
    parser.add_option("-l", "--each-file-line", dest="line",
        default=1000, type='int', 
        metavar="LINE", help="Number line per output file")
    parser.add_option("-n", "--number", dest="number",
        type="int", help="Number of spliited files. This has high priority then <-l>.")
    parser.add_option("-o", "--output-prefix", dest="op",
        metavar="FILEIN", help="Output prefix")
    #parser.add_option("-c", "--choice", dest="choice",
    #    type="choice", choices=["a", "b", "c"], 
    #    default="a", help="Supply an int number")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    lineno = options.line
    number = options.number
    op   = options.op
    if not op:
        op = file 
    op = op.rstrip('.') + '.'
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if number:
        count = ['grep', '-c', '">"', file]
        #print >>sys.stderr, count
        totalsequence = int(check_output(' '.join(count), shell=True))
        #print >>sys.stderr, totalsequence
        lineno = int(totalsequence / number) + 1
        #print >>sys.stderr, lineno
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    count_file = 0
    count_line = 1
    fh_split = ''
    for line in fh:
        if line[0] == '>':
            #print >>sys.stderr, "count_line", count_line
            #print >>sys.stderr, "lineno", lineno
            #print >>sys.stderr, "count_line % lineno", count_line % lineno
            if count_line % lineno == 1:
                #print >>sys.stderr, "count_line", count_line
                count_file += 1
                #print >>sys.stderr, 'here'
                if fh_split:
                    fh_split.close()
                fh_split = open(op+str(count_file), 'w')
                #print >>sys.stderr, fh_split
            count_line += 1
            #print >>sys.stderr, fh_split
            print >>fh_split, line.strip()
        #-------------------------------
        else:
            print >>fh_split, line.strip()
    fh_split.close()
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

if __name__ == '__main__':
    main()
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


