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
    the parameter for Blat. Since the shortest query size that will
    guarantee a match is 2 * stepSize + tileSize -1, the default value
    for this is 2 * 11 + 11 -1 = 32 for nucleotides. Try -stepSize=5
    -tileSize=6 -repMatch=10000 -minScore=0 -fine -minIdentity=0


First blat splits the reference sequence up into " tiles" . The manner
in which it is split depends on two parameters,  -tileSize and
-stepSize,  where -tileSize is the size of the tile and -stepSize
specifies when to start the next tile. The default setting of both for
DNA sequences is 11,  which also means the tiles do not overlap.

For blat to report an alignment,  your query sequence must match at
least two tiles (set via -minMatch) with no mismatches (you can allow
up to one mismatch in the tile by using -oneOff). So if you're trying
to align a 21 bp sequence to a reference using the default setting,
blat will never report an alignment.

To illustrate,  imagine this reference sequence (44bp):

>database
AAAAAAAAAAACCCCCCCCCCCGGGGGGGGGGGTTTTTTTTTTT

and this query sequence (12bp)

>test
GGGGGGGGGGGT:
    
The only way an alignment will be reported is if the tileSize is
set to 1 and the minScore set to less than 12.


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
fasta format or other blat supported formats. Like \
~/home/server-project/2bit/mm9.2bit or \
~/home/server-project/bwa_index/mouse/mm9/mm9.fa")
    parser.add_option("-p", "--parameter-blat", dest="blat_p",
        default="-t=dna -q=rna -minIdentity=100", 
        help="Parameters for blat. Full list of accepted parameter \
would be shown when you hit <blat> in terminal. \
Default <-t=dna -q=rna -minIdentity=100> means target sequence is DNA,\
query sequence is RNA, nomismatch allowed. \
For short sequence (20-30 nt), parameter \
<-t=dna -q=rna -minIdentity=0 -tileSize=6 -stepSize=5 \
-minScore=0 -fine.")
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
    #---Blat alignment-----------
    cmd = ' '.join(['blat', genome, file, output, blat_p])
    os.system(cmd)
    #--Filter mapped results--------
    filter = output + '.filter'
    cmd = ' '.join(['pslCDnaFilter -minId=1 -minCover=1',  
        output, filter])
    os.system(cmd)
    aDict = transferPSLtoBed(open(filter), 0)
    fh = open(filter+'.bed', 'w')
    for valueL in aDict.values():
        print >>fh, '\n'.join(valueL)
    fh.close()
    #-----------end close fh-----------
    inputSet = set()
    for line in open(file):
        if line[0] == '>':
            inputSet.add(line[1:-1])
    mappedSet = set()
    head = 5
    for line in open(output):
        if head:
            head -= 1
            continue
        mappedSet.add(line.split('\t')[9])
    #------------------------------------
    finalSet = set(aDict.keys())
    print >>sys.stderr, "%d FASTA sequences are given in \
file %s, %d sequences can be mapped to given database and \
%d sequences in final output." % \
(len(inputSet), file, len(mappedSet), len(finalSet))
    if inputSet.difference(mappedSet):
        print >>sys.stderr, "The following sequences have \
    not been mapped to database:", '\t'.join(inputSet.difference(mappedSet))
    if mappedSet.difference(finalSet):
        print >>sys.stderr, "The following sequences have \
no perfect match to database, you may want to \
extend their sequence.\n", '\t'.join(mappedSet.difference(finalSet))
    #------------------------------------
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



