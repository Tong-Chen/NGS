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
    this is designed to generate input file for psRobot_mir.

<smRNA seq><Tab><reads1><Tab><reads2>...

'''

import sys
import os
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
        metavar="FILEIN", help="One or multiple FASTA files separated by comma.")
    parser.add_option("-l", "--label", dest="label",
        metavar="FILEIN", help="One or multiple FASTA file names separated by comma. Order must be the same as <-i>.")
    #parser.add_option("-n", "--number", dest="number",
    #    type="int", help="Supply an int number")
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
    fileL = [i.strip() for i in re.split(r'[, ]', file.strip())]
    label = options.label
    labelL = [i.strip() for i in re.split(r'[, ]', label.strip())]
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    aDict = {}
    for file, label in zip(fileL, labelL):
        for line in open(file):
            if line[0] != '>':
                seq = line.strip().strip('N')
                if seq not in aDict:
                    aDict[seq] = {}
                aDict[seq][label] = aDict[seq].get(label, 0) + 1
        #-------------------------------------------------
    #-------------END reading file----------
    for seq, labelD in aDict.items():
        tmpL = [str(labelD.get(label, 0)) for label in labelL]
        print "{}\t{}".format(seq, '\t'.join(tmpL))
    #----close file handle for files-----
    #-----------end close fh-----------
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
    if verbose:
        #print("--Successful %s" % strftime(timeformat, localtime()), file=sys.stderr)
        print >>sys.stderr, "--Successful %s" % strftime(timeformat, localtime())

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    with open('python.log', 'a') as fh:
        #print ("%s\n\tRun time : %s - %s " % \
        #(' '.join(sys.argv), startTime, endTime), file=fh)
        print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


