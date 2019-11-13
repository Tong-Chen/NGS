#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to transfer genome bin count matrix to melted format for plooting.

Input file 

--- data matrix
	ID  A   B   C   D
	chr1_1  1   1   1   1
	chr1_2  1   1   1   1
	chr2_1  1   1   1   1
	chr3_1  1   1   1   1

--- chrom (only first column will be used, so single column file will also OK)
    chr1    0       248956422       chr1
    chr2    0       242193529       chr2
    chr3    0       198295559       chr3
    chr4    0       190214555       chr4
    chr5    0       181538259       chr5
    chr6    0       170805979       chr6
    chr7    0       159345973       chr7
    chr8    0       145138636       chr8

Output (assume bin=100)
value   sample    sample  set
1   1   A   chr1
1   1   B   chr1
1   1   C   chr1
1   1   D   chr1
1   101   A   chr1
1   101   B   chr1
1   101   C   chr1
1   101   D   chr1

Hints: The program will generate a OUT_PFX.summary.xls file containing the figures plotted and their label.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import re
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

def fprint(content):
    """ 
    This is a Google style docs.

    Args:
        param1(str): this is the first param
        param2(int, optional): this is a second param
            
    Returns:
        bool: This is a description of what is returned
            
    Raises:
        KeyError: raises an exception))
    """
    print json_dumps(content,indent=1)

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "\n\t%prog -i file -c chrom -b 50 -o op_prefix -s\n\tzcat file.gz | %prog -i - -c chrom -b 50 -o op_prefix -s"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="Data matrix as described above.")
    parser.add_option("-b", "--bin", dest="bin",
        type="int", 
        help="Bin value (region size) given to `bedtools makewindows -w`.")
    parser.add_option("-c", "--chrom", dest="chrom",
        help="Chromosomes wanted to be analyzed. Chromosome name should be existed in the first column. If unsupplied, all chromosomes including random and unposiitoned ones will be analyzed.")
    parser.add_option("-l", "--labelL", dest="label",
        help="[,] or [ ] separated list of strings to represent sample names..")
    parser.add_option("-p", "--plot-only", dest="plot_only",
        default=False, action="store_true", 
        help="Only do the plot.")
    parser.add_option("-o", "--output-prefix", dest="out_pfx",
        help="Prefix of output files.")
    parser.add_option("-s", "--split", dest="split",
        default=False, action="store_true", help="Output each chromosome into different files. Default False.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def plot(file, chr, label_legend):
    if not os.path.exists(file):
        return
    if chr == 'all':
        cmd = ['s-plot scatterplot2 -X pos -Y value -c sample -S sample', 
               '-x "Chromosome positions" -y "Coverage" -A 0.5', 
               '-f', file, '-I', label_legend, 
               '-F \'+facet_wrap(~chr, ncol=1, scale="free")\'', 
               '-w 24 -u 12 -E png']
    else:
        cmd = ['s-plot scatterplot2 -X pos -Y value -c sample -S sample', 
               '-x "Chromosome positions" -y "Coverage" -A 0.5', 
               '-f', file, '-I', label_legend, 
               '-w 24 -u 12 -E png']
    if debug:
        print ' '.join(cmd)
    os.system(' '.join(cmd))

#-------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    bin  = options.bin
    chrom = options.chrom
    if chrom:
        chromL = [line.split()[0] for line in open(chrom)]
    else:
        chromL = []
    out_pfx = options.out_pfx
    plot_only = options.plot_only
    split = options.split
    verbose = options.verbose
    labelL = re.split(r'[, ]*', options.label)
    labelL = '"'+','.join(["'"+i+"'" for i in labelL])+'"'
    global debug
    debug = options.debug
    #-----------------------------------
    if plot_only:
        if options.split:
            for chr in chromL:
                file = out_pfx + '.'+chr+'.melt.xls'
                plot(file, chr, labelL)
        else:
            out = out_pfx +'.melt.xls'
            plot(out, 'all', labelL)
        return 
    #--------------------------------------------        
    summary = out_pfx + '.summary.xls'
    summary_fh = open(summary, 'w')
    if not options.split:
        out = out_pfx +'.melt.xls'
        out_fh = open(out, 'w')
        print >>out_fh, "value\tpos\tsample\tchr"
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    '''
    filefhD = {'chr1': [file_name, file_fh], 'chr2':[file_name, file_fh]}
    '''
    filefhD = {}
    header = 1
    for line in fh:
        lineL = line.split()
        if header:
            headerL = lineL[1:]
            header -= 1
            continue
        ID = lineL[0]
        chr, pos = ID.rsplit('_', 1)
        if chromL and chr not in chromL:
            continue
        coord = str((int(pos)-1)*bin+1)
        valueL = lineL[1:]
        tmpL = []
        if not options.split:
            for value, sample in zip(valueL, headerL):
                tmpL.append('\t'.join([value, coord, sample, chr]))
            print >>out_fh, '\n'.join(tmpL)
        else:
            for value, sample in zip(valueL, headerL):
                tmpL.append('\t'.join([value, coord, sample]))
            if chr not in filefhD:
                file = out_pfx + '.'+chr+'.melt.xls'
                fh_tmp = open(file, 'w')
                filefhD[chr] = [file, fh_tmp]
                print >>fh_tmp, "value\tpos\tsample"
            fh_tmp = filefhD[chr][1]
            print >>fh_tmp, '\n'.join(tmpL)
    #-------------END reading file----------
    if not options.split:
        out_fh.close()
    else:
        for file, fh_tmp in filefhD.values():
            fh_tmp.close()

    if not options.split:
        plot(out, 'all', labelL)
        print >>summary_fh, "\t".join([out+'.scatterplot2.png', 'all'])
    else:
        for chr, fhL in filefhD.items():
            file = fhL[0]
            plot(file, chr, labelL)
            print >>summary_fh, "\t".join([file+'.scatterplot2.png', chr])
    #----close file handle for files-----
    summary_fh.close()
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


