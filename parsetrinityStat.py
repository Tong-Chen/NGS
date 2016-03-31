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
    This is used to transfer the output of TrinityStats.pl to pandoc 
    markdown table.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

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
        metavar="FILEIN", help="The output of TrinityStats.pl")
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
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    line = fh.readline()
    while line[0] != '#': #Get the first # line
        line = fh.readline()
    line = fh.readline() #Get the ## line as title
    title = line.replace('#', '').strip()
    line = fh.readline() # Skip the line with all #
    line = fh.readline() # Get the next line
    print 'Table:', title
    print 
    print '|'.join(['Item', 'Count'])
    print ':--------------------------|------------------------:'
    while line[0] != '#':  # Stop when meets the next line with all '#'
        if line != '\n':
            key, value = [i.strip() for i in line.split(':')]
            print '|'.join([key, value])
        line = fh.readline()
    #---------------------------------------
    #----END summary--------------------
    print 
    print 'Table:', title
    print 
    print '|'.join(['Item', 'Count'])
    print ':--------------------------|------------------------:'
    line = fh.readline()
    title = line.strip().replace('#', '').strip(' :')
    line = fh.readline()
    line = fh.readline()
    while line[0] != '#':
        if line != '\n':
            key, value = [i.strip() for i in line.split(':')]
            print '|'.join([key, value])
        line = fh.readline()
    print ':', title
    #----END summary--------------------

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


