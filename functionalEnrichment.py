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
    This is designed to do functional enrichment for given entrez ids or
    groups od entrez IDs.

Normally input file is:
    839418  cbf0._vs_.cbf3_up
    837073  cbf0._vs_.cbf3_up
    837379  cbf0._vs_.cbf3_up
    837552  cbf0._vs_.cbf3_up
    837607  cbf0._vs_.cbf3_up
    837714  cbf0._vs_.cbf3_up
    2745749 cbf0._vs_.cbf3_dw
    837980  cbf0._vs_.cbf3_dw
    838005  cbf0._vs_.cbf3_dw
    838152  cbf0._vs_.cbf3_dw
    838153  cbf0._vs_.cbf4_up
    838554  cbf0._vs_.cbf4_up
    838870  cbf0._vs_.cbf4_dw
    838947  cbf0._vs_.cbf4_dw
    839223  cbf0._vs_.cbf4_dw

or 
    839418
    837073
    837379
    837552
    837607
    837714


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
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A multiple column file as \
described above.")
    parser.add_option("-s", "--group_col", dest="group_col",
        default=2, help="Specify the columns used to group. \
Default 2 indicating group files based on the content of \
the second column.")
    parser.add_option("-o", "--output_col", dest="output_col",
        default=1, help="Specify the columns need to be output. \
Default 1 indicating output the first column.")
    parser.add_option("-g", "--genome_label", dest="genome",
        help="Species name for geneontology analysis")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def output_func(aDict, file):
    fileD = {}
    for key, valueL in aDict.items():
        out_file = file + '.' + key
        fileD[key] = out_file
        out_fh = open(out_file, 'w')
        print >>out_fh, '\n'.join(valueL)
        out_fh.close()
    return fileD
#----END of output-------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    group_col = int(options.group_col) - 1
    output_col = int(options.output_col) - 1
    genome = options.genome
    output_dir = os.path.dirname(file)
    #os.system("mkdir -p "+output_dir)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    for line in fh:
        lineL = line.strip().split('\t')
        output = lineL[output_col]
        group = lineL[group_col]
        if group not in aDict:
            aDict[group] = []
        else:
            aDict[group].append(output)
    #-------------END reading file----------
    fileD = output_func(aDict, file)
    #-----Perform GO analysis-----------------
    for key, file in fileD.items():
        print >>sys.stderr, "Performing functional enrichment for",file
        cmd = ["findGO.m.pl", file, genome, output_dir+'/'+key, "-cpu 20"] 
        print >>sys.stderr,  ' '.join(cmd)
        os.system(' '.join(cmd))
        cat_cmd = ['cat', output_dir+'/'+key+"/*.xls", "| cut -f 1-4,6",
                "| grep -v", '^TermID', 
                '| sed "1 i\TermID\tTerm\tFDR\tneg_log10FDR\tcount"', 
                '>', output_dir+'/'+key+'/mergeForWC']
        print >>sys.stderr,  ' '.join(cat_cmd)
        os.system(' '.join(cat_cmd))
        wc = ["wordcloud.sh -f", output_dir+'/'+key+'/mergeForWC', 
                "-w Term -F neg_log10FDR"]
        print >>sys.stderr, ' '.join(wc)
        os.system(' '.join(wc))
    os.system("/bin/rm -f *.fg.tmp")
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


