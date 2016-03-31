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
    This is designed to add extra columns with specified values to
    separate several groups of data.
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
        metavar="FILEIN", help="A multiple column file")
    parser.add_option("-n", "--num-data-col", dest="num_data_col",
        help="The number of columns for each group.")
    parser.add_option("-N", "--num-label-col", dest="num_label_col",
        help="The number of label columns before data columns.")
    parser.add_option("-e", "--num-extra-col", dest="num_extra_col",
        default=10, help="The number of extra columns to add. Default 10.")
    parser.add_option("-v", "--val-extra-col", dest="val_extra_col",
        default='0', help="The values in extra columns.")
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
    num_data_col = int(options.num_data_col)
    num_label_col = int(options.num_label_col)
    num_extra_col = int(options.num_extra_col)
    val_extra_col = options.val_extra_col
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    header = 1
    for line in fh:
        lineL = line.split('\t')
        if header:
            lenLine = len(lineL)
            dataI = []
            lastIndex = lenLine - num_data_col
            for index in range(num_label_col, lastIndex, num_data_col):
                dataI.append((index, index+num_data_col)) 
            header -= 1
            extra_col = '\t'.join([val_extra_col] * num_extra_col)
            lastIndex = lenLine - num_data_col
        #-------------------------------------
        dataLineL = ['\t'.join(lineL[index[0]:index[1]])+'\t'+extra_col \
            for index in dataI]
        print "%s\t%s\t%s" % ('\t'.join(lineL[:num_label_col]),
                '\t'.join(dataLineL), '\t'.join(lineL[lastIndex:])) 
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


