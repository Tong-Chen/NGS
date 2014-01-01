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
    This is designed to trim Sequences in FASTA file.
'''

import sys
import os
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
        metavar="FILEIN", help="A FATSA file.")
    parser.add_option("-s", "--start", dest="start",
        metavar="0", default=0, help="The position of the first \
letter kept. \
Default 0 indicates no trimming of left part. All positive \
numebr is OK.")
    parser.add_option("-e", "--end", dest="end",
        metavar="-1", default=-1, help="The position of the \
last letter kept. Default -1 indicates trimming the last letter. \
Any number is accepted.")
    parser.add_option("-l", "--min-len", dest="min_len",
        metavar="16", default=16, help="Any reads shorted \
than this value will be discared.")
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
    file  = options.filein
    start = int(options.start)
    end   = int(options.end)
    min_len = int(options.min_len)
    if end > 0:
        assert end > start, "End smaller than start"
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    key = ''
    for line in fh:
        if line[0] == '>':
            if key:
               seq = ''.join(seqL)
               trim_seq = seq[start:end]
               if len(trim_seq) >= min_len:
                   print "%s%s" % (key, trim_seq)
            key = line
            seqL = []
        else:
           seqL.append(line.strip())
    #-------------END reading file----------
    #----The last one if exists---------
    if key:
       seq = ''.join(seqL)
       trim_seq = seq[start:end]
       if len(trim_seq) >= min_len:
           print "%s%s" % (key, trim_seq)
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



