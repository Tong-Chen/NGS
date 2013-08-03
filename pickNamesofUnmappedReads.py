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

desc='''
Functionla description

This is designed to get the names of unmapped reads if its mate reads
mapped.


Output format: (1 and 2 represents two ends, 
          flag 69 given 1, flag 133 given 2)
reads_name 1
reads_name 2
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
        metavar="FILEIN", help="A sam file or a pipeline from \
samtools view")
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
    flagD = {'69':'1', '133':'2'}
    aDict = {}
    for line in fh:
        lineL = line.split('\t',2)
        name = lineL[0]
        flag = lineL[1]
        if flag == '69' or flag == '133':
            flag = flagD[flag]
        if name not in aDict:
            aDict[name] = [flag]
        else:
            aDict[name].append(flag)
        #-------------------------------------
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    for key, valueL in aDict.items():
        if len(valueL) == 1:
            print "%s\t%s" % (key, valueL[0])
        elif len(valueL) == 2:
            #pass
            print "%s\t%s" % (key, valueL[0])
            print "%s\t%s" % (key, valueL[1])
        else:
            print >>sys.stderr, "Unexpected things %s" \
                % '\t'.join(valueL)
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



