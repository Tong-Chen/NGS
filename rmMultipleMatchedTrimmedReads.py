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
'''
Functionla description

This is designed to remove reads from tophat output 
which have multiple matched locations
caused by base trimming.

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
    desc = ""
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A file containing \
the names of trimmed reads and the ends they belong to. \
Such as reads_name1\t1\nreads_name2\t2")
    parser.add_option("-s", "--sam", dest="sam",
        metavar="sam", help="Output of samtools view -X \
accepted_hits.bam, accept - which means sys.stdin")
    parser.add_option("-m", "--max-num-hits", dest="mnh",
        default=1, metavar=1, help="The number of allowed multiple \
mapped regions. Default 1")
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
    sam = options.sam
    mnh = int(options.mnh)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    leftD = {}
    rightD = {}
    for line in open(file):
        name, flag = line.split()
        if flag == '1':
            leftD[name] = flag
        elif flag == '2':
            rightD[name] = flag
    #-------------END reading file----------
    if sam == '-':
        fh_sam = sys.stdin
    else:
        fh_sam = open(sam)
    #--------------------------------
    for line in fh_sam:
        output = 1
        lineL = line.split()
        key = lineL[0]
        flag = lineL[1]
        if (flag.find('1') != -1 and key in leftD) \
            or (flag.find('2') != -1 and key in rightD):
            for i in lineL[11:]:
                if i.find('NH:i:') == 0:
                    hit = int(i[5:])
                    if hit > mnh:
                        output = 0
                    break
            #----------------------------
        #-----------------------------------------
        if output:
            print line,
    #----close file handle for files-----
    if sam != '-':
        fh_sam.close()
    #-----------end close fh-----------
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



