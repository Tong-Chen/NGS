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
    This is designed to get the bed position of gene or mRNA
    from Cufflink GTF file. Only exon type records will be used.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from bs4 import BeautifulSoup

#reload(sys)
#sys.setdefaultencoding('utf8')

#from multiprocessing.dummy import Pool as ThreadPool

debug = 0

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
        metavar="FILEIN", help="GTF file with at least \
<exon> record.")
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
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    posD = {}
    for line in fh:
        lineL = line.split("\t")
        type = lineL[2]
        if type != "exon":
            continue
        start = int(lineL[3]) - 1
        end   = int(lineL[4])
        #print start,end
        strand = lineL[6]
        chr = lineL[0]
        nameD = dict([i.split() \
            for i in lineL[8].replace('"','').split("; ")])
        gene_id = nameD['gene_id']
        transcript_id = nameD['transcript_id']
        #print gene_id, transcript_id
        if gene_id not in posD:
            posD[gene_id] = [chr, start, end, gene_id, '.', strand]
        else:
            if start < posD[gene_id][1]:
                posD[gene_id][1] = start
            if end > posD[gene_id][2]: 
                posD[gene_id][2] = end
        #----------------------------------------------------
        if transcript_id not in posD:
            posD[transcript_id] = [chr, start, end, transcript_id, '.', strand]
        else:
            if start < posD[transcript_id][1]:
                posD[transcript_id][1] = start
            if end > posD[transcript_id][2]: 
                posD[transcript_id][2] = end
        #----------------------------------------------------
    #-------------END reading file----------
    for value in posD.values():
        value[1] = str(value[1])
        value[2] = str(value[2])
        print '\t'.join(value)
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


