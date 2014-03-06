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

'''

import sys
import os
import random
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

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
        metavar="FILEIN", help="Output of compareMotifs.pl")
    parser.add_option("-n", "--number-motif", dest="number",
        metavar="NUM", help="The number of motifs you generated \
for each simulation. Since generated motifs will be done similarity \
check between each other, you may not get so many motifs sometimes.")
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
    num = int(options.number)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    group = 0
    skip_additional_header = 1
    for line in fh:
        if line.find('Determining similar motifs') != -1 \
            and line.find(' reduced to ') != -1:
            #print line
            skip_additional_header = 0
            if group == 0:
                print "group\treal\tcount_test\tselected_count\tcount_match"
            if group and real == count_test:
                print "%d\t%d\t%d\t%d\t%d\t%f" % \
                    (group, real, count_test, selected_count,
                    count_match, count_match*1.0/selected_count)
            group += 1
            #print real, count_test
            real = int(line.strip().split()[-2])
            if verbose:
                print " ###", real 
            if real >= num:
                process = 1
                random_time = real - num
            else:
                process = 0
                random_time = 0
            count_test = 0
            selected_count = 0
            count_match = 0
            skip = 0 
            if verbose:
                print real, num, random_time
        elif skip_additional_header == 0 and process == 1 \
                and line.find(' similar to') != -1:
            count_test += 1
            if skip == 0:
                select = 1
                if random_time != 0:
                    select = random.randint(0, 1)
                    if select == 0:
                        random_time -= 1
                if select:
                    selected_count += 1
                    if not line.rstrip().endswith("similar to"):
                        count_match += 1
                if selected_count == num:
                    #skip additional ones
                    skip = 1
            #------------------------------------
            if verbose:
                print line,
                print select, random_time, real, count_test, selected_count, count_match
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
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



