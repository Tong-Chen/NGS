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

This program is designed to locate a mRNA fragment to genome. It first
maps mRNA fragment to genome using blat and then transfers the output
to bed file.

Potential bug:
    1. If a fragment or one of the spliced fragments is too short, it can not
    be mapped to genome with default parameter. You may want to modify
    the parameter for Blat.
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
from transferPSLtoBed import transferPSLtoBed

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
        metavar="FILEIN", help="A fatsa file contains the mRNA \
fragments.")
    parser.add_option("-g", "--genome-fasta", dest="genome",
        metavar="GENOME", help="A genome assembl file in \
fasta format or other blat supported formats.")
    parser.add_option("-p", "--parameter-blat", dest="blat_p",
        default="-t=dna -q=rna -minIdentity=100", 
        help="Parameters for blat. Full list of accepted parameter \
would be shown when you hit <blat> in terminal. \
Default <-t=dna -q=rna -minIdentity=100> means target sequence is DNA,\
query sequence is RNA, nomismatch allowed.")
    parser.add_option("-o", "--output-prefix", dest="op",
        help="The prefix for output file; \
Default the string given to -i.")
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
    genome = options.genome
    blat_p = options.blat_p
    output = options.op
    if not output:
        output = file
    output = output + '.blat.psl'
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    cmd = ' '.join(['blat', genome, file, output, blat_p])
    os.system(cmd)
    aDict = transferPSLtoBed(open(output), 5)
    fh = open(output+'.bed', 'w')
    for valueL in aDict.values():
        print >>fh, '\n'.join(valueL)
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



