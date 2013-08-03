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
Functionla description

This is used to trim pair-end reads in fastq format. 

Those reads which name in the file given by -i will be trimmed. If you
want to trim all reads,  an empty file cound be given.
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    global desc
    print >>sys.stderr, desc
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A file containing \
the names of reads you want to trim and the ends they belong to. \
Such as reads_name1\t1\nreads_name2\t2")
    parser.add_option("-p", "--prefix", dest="prefix",
        metavar="reads", help="For file name <reads_1.fq>, <reads> \
would be suitable here. <_1> and <_2> will be added automatically.")
    parser.add_option("-s", "--suffix", dest="suffix",
        default="fq", help="Usually <fq> or <fastq>")
    parser.add_option("-f", "--first", dest="first_base",
        default=1, help="First base to keep (1-based). Default is 1.")
    parser.add_option("-l", "--last", dest="last_base",
        default="", help="Last base to keep (1-based). Default entire end.")
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
    prefix = options.prefix
    suffix = options.suffix
    left_fq = options.prefix + '_1.' + suffix
    right_fq = options.prefix + '_2.' + suffix
    left_out_fq = options.prefix + '.trim_1.' + suffix
    right_out_fq = options.prefix + '.trim_2.' + suffix
    first_base = int(options.first_base) - 1
    if options.last_base:
        last_base = int(options.last_base) + 1
    else:
        last_base = options.last_base
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    leftD = {}
    rightD = {}
    for line in fh:
        name, flag = line.split()
        if flag == '1':
            leftD[name] = flag
        elif flag == '2':
            rightD[name] = flag
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    fh1 = open(left_fq)
    fh1_out = open(left_out_fq, 'w')
    left_i = 0
    for line in fh1:
        left_i += 1
        mod = left_i % 4
        if mod == 1:
            print >>fh1_out, line,
            start = 0
            end = '' # A random value
            name = line.split()[0][1:]
            if name in leftD:
                start = first_base
                end = last_base
        elif mod == 2 or mod == 0:
            if end != '':
                print >>fh1_out, line.rstrip()[start:end]
            else:
                print >>fh1_out, line.rstrip()[start:]
        else:
            print >>fh1_out, line,
        #------------------------------------
    #----------------------------------------
    fh1.close()
    fh1_out.close()
    fh2 = open(right_fq)
    fh2_out = open(right_out_fq, 'w')
    #------------------------------------------
    right_i = 0
    for line in fh2:
        right_i += 1
        mod = right_i % 4
        if mod == 1:
            print >>fh2_out, line,
            start = 0
            end = '' # A random value
            name = line.split()[0][1:]
            if name in rightD:
                start = first_base
                end = last_base
        elif mod == 2 or mod == 0:
            if end != '':
                print >>fh2_out, line.rstrip()[start:end]
            else:
                print >>fh2_out, line.rstrip()[start:]
        else:
            print >>fh2_out, line,
        #------------------------------------
    #----------------------------------------
    fh2.close()
    fh2_out.close()
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



