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
    This is designed to find the reason why reads are unmapped and
    what are they mapping to.
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
    parser.add_option("-f", "--fasta", dest="fasta",
        help="The fasta file used for blast analysis.")
    parser.add_option("-n", "--num_seq", dest="num_seq",
        type="int", help="Number of fasta sequences used for downstream analysis. \
Ignore this parameter if all sequences are used.")
    parser.add_option("-l", "--label", dest="label",
        help="The label for this analysis.")
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="The output of blastn. Normally \
the format is \
<6 qseqid sseqid qlen qstart qend evalue bit score \
pident nident sstrand sscinames>")
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
    label = options.label
    num_seq = options.num_seq
    fasta = options.fasta
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    total_seq = 0
    for line in open(fasta):
        if line[0] == '>':
            total_seq += 1
    #-----------------------------------
    if num_seq and total_seq > num_seq:
        total_seq = num_seq
    mapped = 0
    staD = {}
    map_percent_threshold = 0.4
    map_ident_threshold = 0.8
    #---------------
    key = ''
    outputL = []
    for line in fh:
        lineL = line.strip().split('\t')
        newkey = lineL[0]
        if key and key != newkey:
            if outputL[5] >= map_ident_threshold and match_percent >= map_percent_threshold:
                staD[outputL[-1]] = staD.get(outputL[-1], 0)+1
            #outputL = [str(i) for i in outputL]
            #print '\t'.join(outputL)
            outputL = []
        key = newkey
        sid = lineL[1]
        qlen = int(lineL[2])
        qstart = int(lineL[3])
        #qend   = int(lineL[4]) - qlen - 1
        qend = int(lineL[4])
        match_percent = (qend-qstart+1)*1.0/qlen
        pident  = float(lineL[7])
        nident  = int(lineL[8])
        sscinames = lineL[-1] 
        #total_pident = nident * 1.0 / qlen * 100
        if (len(outputL)==0) or (outputL[5]<pident):
            outputL = [key, sid, qlen, qstart, qend, pident, match_percent, sscinames]
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    mapped = sum(staD.values())
    unmapped = total_seq - mapped
    print "\t".join(["Unmapped", str(unmapped), label])
    keyL = staD.keys()
    keyL.sort(key=lambda x: staD[x], reverse=True)
    for key in keyL[:10]:
        print "\t".join([str(i) for i in [key,staD[key],label]])
        staD.pop(key)
    other = sum(staD.values())
    print "\t".join(["Other_species", str(other), label])

    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
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
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


