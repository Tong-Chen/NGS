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
    This is designed to generate random sequences bu given alphabets.
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import random

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-alphabet", dest="alpha",
        metavar="ALPHA", default="ACGTRYMKWSBDHVN", help="\
Accept a string containing all alphabets one needed to \
construct random string. Default 'ACGTRYMKWSBDHVN'.")
    parser.add_option("-p", "--prefix", dest="prefix",
        help="The prefix for each FASTA sequence in output gile. \
Normally the sequence name will be in the format <prefix>_<num>.")
    parser.add_option("-n", "--num-of-seq", dest="count",
        help="The numebr of random sequence you want to generate.")
    parser.add_option("-l", "--len-of-seq", dest="len_seq",
        help="The length of random sequence you want to generate. \
Fixed length string would be generated when single number is given. \
Random length string with length varies form given range will be \
generated when two numbers supplied with ',' like '5,10'.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.count != None, "A count number needed for -n"
    assert options.len_seq != None, "A len_seq number needed for -l"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    alpha = options.alpha
    prefix = options.prefix
    count = int(options.count)
    len_seq = [int(i) for i in options.len_seq.split(',')]
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    for i in range(count):
        print ">%s_%s" % (prefix, str(i))
        random_len = random.randint(len_seq[0],len_seq[-1])
        print ''.join(random.sample(alpha, random_len))
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



