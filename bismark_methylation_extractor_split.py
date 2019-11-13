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
    This is designed to split bismark_methylation_extractor output file to be suitable for following analysis.

Two types of input file are allowed.

The first type is the output of bismark_methylation_extractor or other format specified below (without the header line):

#<chromosome> <position> <strand> <count methylated> <count non-methylated> <C-context> <trinucleotide context>

The  second type would be in format specified below (with header line) 
#Chr    Pos     Ref     Chain   Total   Met     UnMet   MetRate Ref_context     Type
chr1    14485   C       +       1       1       0       1       CCG     CHG
chr1    14486   C       +       1       1       0       1       CGT     CpG
chr1    14489   C       +       1       1       0       1       CCC     CHH

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import gzip
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
    usages = "\n\t%prog -i file\n\tzcat gzippedfile | %prog -i -"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", 
        help="One of the two types of files in format as described above. Gzipped file should be input as STDIN (-).")
    parser.add_option("-l", "--label", dest="label",
        help="One simple string normally sample name as label.")
    parser.add_option("-t", "--type", dest="type",
        default='first', help="<first> or <second>. Default <first>.")
    parser.add_option("-s", "--statistics-only", dest="statistics",
        default=False, help="Compute statistics, like distribution of reads coverage, methylation rate. Default <False>.")
    parser.add_option("-o", "--output-prefix", dest="op_prefix",
        help="Prefix of output files.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def bismark_split(fh, strandD, fhD, countD, methylRateD, statistics_only):
    #<chromosome> <position> <strand> <count methylated> <count non-methylated> <C-context> <trinucleotide context>
    for line in fh:
        lineL = line.strip().split()
        methyl = int(lineL[3])
        unmethyl = int(lineL[4])
        total = methyl+unmethyl
        if total == 0:
            continue
        type = lineL[5]
        if type == "CG":
            type = "CpG"
        countD[type] += 1
        freqC = "{:.2f}".format(100.0*methyl/total)
        if debug:
            print >>sys.stderr, "\t".join([str(methyl), str(total), freqC, type])
        methylRateD[type][freqC] = methylRateD[type].get(freqC, 0)+1
        if not statistics_only:
            chr = lineL[0]
            base = lineL[1]
            strand = strandD[lineL[2]]
            freqT = "{:.2f}".format(100.0*unmethyl/total)
            chrBase = '.'.join([chr, base, lineL[6]])
            print >>fhD[type], '\t'.join([chrBase, chr, base, strand, str(total), freqC, freqT])
    #-------------------------------------
#----------------------------------------------
def second_split(fh, strandD, fhD, countD, methylRateD, statistics_only):
    #Chr    Pos     Ref     Chain   Total   Met     UnMet   MetRate Ref_context     Type
    header = 1
    for line in fh:
        if header:
            header -= 1
            continue
        lineL = line.strip().split()
        total = lineL[4]
        if total == '0':
            continue
        type = lineL[-1]
        if type == "CG":
            type = "CpG"
        countD[type] += 1
        MetRate = float(lineL[7])*100
        freqC = "{:.2f}".format(MetRate)
        methylRateD[type][freqC] = methylRateD[type].get(freqC, 0)+1
        if not statistics_only:
            strand = strandD[lineL[3]]
            chr = lineL[0]
            base = lineL[1]
            freqT = "{:.2f}".format(100-MetRate)
            chrBase = '.'.join([chr, base, lineL[8]])
            print >>fhD[type], '\t'.join([chrBase, chr, base, strand, str(total), freqC, freqT])
    #-------------------------------------
#----------------------------------------------
def output(countD, methylRateD, label, op_prefix):
    '''
    countD = {"CpG":0, "CHH":0, "CHG":0}
    methylRateD = {"CpG":{}, "CHH":{}, "CHG":{}}
    '''
    count_fh = open(op_prefix+'.c_count.xls', 'w')
    methyl_fh = open(op_prefix+'.methylRate.xls', 'w')
    print >>count_fh, "value\tvariable\tsample"
    for key, value in countD.items():
        print >>count_fh, "\t".join([str(value), key, label])
    count_fh.close()
    
    print >>methyl_fh, "value\tvariable\tsample\ttype"
    for type, valueD in methylRateD.items():
        for variable, value in valueD.items():
            print >>methyl_fh, '\t'.join([str(value), variable, label, type])
    methyl_fh.close()
#---------------------------------------------------
def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    label = options.label
    type = options.type
    statistics_only = options.statistics
    op_prefix = options.op_prefix
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------

    strandD = {"+":"F", "-":"R"}
    fhD = {}
    if not statistics_only:
        CpG = op_prefix+'.CpG.txt.gz'
        CHG = op_prefix+'.CHG.txt.gz'
        CHH = op_prefix+'.CHH.txt.gz'
        fhD = {'CpG': gzip.open(CpG, 'wb'), 
               'CHH': gzip.open(CHH, 'wb'), 
               'CHG': gzip.open(CHG, 'wb')} 
        for tmp_fh in fhD.values():
            print >>tmp_fh, "chrBase	chr	base	strand	coverage	freqC	freqT"
    countD = {"CpG":0, "CHH":0, "CHG":0}
    methylRateD = {"CpG":{}, "CHH":{}, "CHG":{}}

    if type == 'second':
        second_split(fh, strandD, fhD, countD, methylRateD, statistics_only)
    elif type == 'first':
        #<chromosome> <position> <strand> <count methylated> <count non-methylated> <C-context> <trinucleotide context>
        bismark_split(fh, strandD, fhD, countD, methylRateD, statistics_only)
    else:
        print >>sys.stderr, "Unsupport type"
        sys.exit(1)
    if not statistics_only:
        for tmp_fh in fhD.values():
            tmp_fh.close()
        #---------------------------------------------
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    output(countD, methylRateD, label, op_prefix)
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


