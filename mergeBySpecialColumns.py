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
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
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
        metavar="FILEIN", help="Files sorted by target columns.")
    parser.add_option("-c", "--columns", dest="column",
        default='1', help="Specify columns used as keys. Default 1 meaning first column. Accept one number or '1,2' format numbers.")
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
    cols = [int(i)-1 for i in options.column.split(',')]
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    header = 1
    aDict = {}
    key = ''
    for line in fh:
        line = line.strip('\n')
        lineL = line.split('\t')
        if header:
            header -= 1
            print line
            len_lineL = len(lineL)
            continue
        #---------------------------------------
        if key and key != lineL[0]:
            print >>sys.stderr, 'Output ' + key
            #print >>sys.stderr, aDict[key]
            print '\t'.join(['|'.join(i) for i in aDict[key]])
            aDict = {}
        key = lineL[0]
        #print >>sys.stderr, key
        tmpL = lineL[:]
        len_tmpL = len(tmpL)
        #if len_tmpL < len_lineL:
        #    print >>sys.stderr, tmpL
        for i in range(len_lineL):
            if i >= len_tmpL:
                tmpL.append(['NA'])
                continue
            if not tmpL[i]:
                tmpL[i] = ['NA']
            else:
                tmpL[i] = [tmpL[i]]
        #--------------------------------------
        if key not in aDict:
            aDict[key] = tmpL
            #print tmpL
        else:
            for i in range(len_lineL):
                old = aDict[key][i] 
                new = tmpL[i][0]
                if new != 'NA' and new not in old:
                    if old == 'NA':
                        aDict[key][i].remove('NA')
                    #aDict[key][i].remove('NA')
                    aDict[key][i].append(new)
                    #----------------
                #----------------
            #print aDict[key]
        #---------------------------------------
    #-------------END reading file----------
    #print key
    if key:
        print >>sys.stderr, 'Output ' + key
        #print >>sys.stderr, aDict[key]
        print '\t'.join(['|'.join(i) for i in aDict[key]])
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


