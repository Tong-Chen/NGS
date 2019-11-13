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
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="Bismark aligned BAM file. REQUIRED.")
    parser.add_option("-L", "--file-label", dest="label",
        metavar="LABEL", help="File label")
    parser.add_option("-S", "--summary-only", dest="summary",
        default=0, type='int', metavar="Upper S", help="Specify 1 to only do summary part.")
    parser.add_option("-o", "--output-prefix", dest="op",
        metavar="FILEIN", help="Prefix of output files. REQUIRED.")
    parser.add_option("-g", "--genome-fa", dest="genome",
        metavar="GENOME FASTA", help="Genome sequence file. REQUIRED.")
    parser.add_option("-f", "--methylC-count-for-all-reads", dest="full",
        default="", metavar="FILEIN", 
        help="Called methylated C sites using all reads. Normally <SingleCmet.gz>.")
    parser.add_option("-l", "--percentile-floor", dest="per_lb",
        metavar="PERCENTILE_LOW_BOUND", default=5, type="int", 
        help="Sampling starts from this percentile. An integer between 0 and 100. default=5")
    parser.add_option("-u", "--percentile-ceiling", dest="per_ub",
        metavar="PERCENTILE_UP_BOUND", default=95, type="int", 
        help="Sampling ends this percentile. An integer between 0 and 100. default=95. ")
    parser.add_option("-s", "--percentile-step", dest="per_step",
        metavar="PERCENTILE_STEP", default=5, type="int", 
        help="Sampling frequency. Smaller value means more sampling times. An integer between 0 and 100. default=5.")
    parser.add_option("-m", "--mapped-reads", dest="map_reads",
        help="Bismark output statistics or a number")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def methylC(output, genome):
    prefix = output[:-4]
    pileup = ['samtools mpileup -f', genome, output, '>', prefix+'.pileup']
    pileup = ' '.join(pileup)
    if debug:
        print >>sys.stderr, pileup
    else:
        os.system(pileup)
    metlevel = ['SingleC_MetLevel.pl', genome, prefix+'.pileup | gzip -c', '>', prefix+'.SingleCmet.gz']
    metlevel = ' '.join(metlevel)
    if debug:
        print >>sys.stderr, metlevel
    else:
        os.system(metlevel)
#--------------------END methylC-------------------

def downSample2(file, per_lb, per_ub, per_step, op, genome=''):
    for per in range(per_lb, per_ub+per_step, per_step):
        per = str(per / 100)
        output = op+'_'+per+'.bam'
        downsample = ["samtools view -b -s", per, file, '>', output]
        downsample = ' '.join(downsample)
        if debug:
            print >>sys.stderr, downsample
        else:
            os.system(downsample)
        methylC(output, genome)

def downSample1(file, per_lb, per_ub, per_step, op, genome=''):
    for per in range(per_lb, per_ub+per_step, per_step):
        per = str(per / 100)
        prefix = op+'_'+per
        if os.path.exists(prefix+'.SingleCmet.gz'):
            print >>sys.stderr, "Exists gz "+prefix
            continue
        elif os.path.exists(prefix+'.SingleCmet'):
            print >>sys.stderr, "Compressing "+prefix 
            os.system("gzip "+prefix+'.SingleCmet')
            continue
        downsample = ["samtools view -u -s", per, file, '|', 
            'samtools mpileup -f', genome, '- >', prefix+'.pileup']
        downsample = ' '.join(downsample)
        if debug:
            print >>sys.stderr, downsample
        else:
            os.system(downsample)
        #methylC(output, genome)
        metlevel = ['SingleC_MetLevel.pl', genome, prefix+'.pileup | gzip -c', '>', prefix+'.SingleCmet.gz']
        metlevel = ' '.join(metlevel)
        if debug:
            print >>sys.stderr, metlevel
        else:
            os.system(metlevel)
            rm = "/bin/rm -f "+prefix+'.pileup &'
            os.system(rm)
#----------------downSample----------------------------
def downSample(file, per_lb, per_ub, per_step, op, genome=''):
    for per in range(per_lb, per_ub+per_step, per_step):
        per = str(per / 100)
        prefix = op+'_'+per
        if os.path.exists(prefix+'.SingleCmet.gz'):
            print >>sys.stderr, "Exists gz "+prefix
            continue
        elif os.path.exists(prefix+'.SingleCmet'):
            print >>sys.stderr, "Compressing "+prefix 
            os.system("gzip "+prefix+'.SingleCmet')
            continue
        downsample = ["samtools view -u -s", per, file, '|', 
            'samtools mpileup -f', genome, '- |', 
            'SingleC_MetLevel.pl', genome, '- | gzip -c', '>', prefix+'.SingleCmet.gz']
        downsample = ' '.join(downsample)
        if debug:
            print >>sys.stderr, downsample
        else:
            os.system(downsample)
        #methylC(output, genome)
        #metlevel = ['SingleC_MetLevel.pl', genome, prefix+'.pileup | gzip -c', '>', prefix+'.SingleCmet.gz']
        #metlevel = ' '.join(metlevel)
        #if debug:
        #    print >>sys.stderr, metlevel
        #else:
        #    os.system(metlevel)
        #    rm = "/bin/rm -f "+prefix+'.pileup &'
        #    os.system(rm)
#----------------downSample----------------------------

def countLines(file, typeL):
    countD = {}
    for type in typeL:
        countD[type] = 0
    with gzip.open(file, 'rb') as f:
        for line in f:
            line = line.rstrip()
            for type in typeL:
                if line.endswith(type):
                    countD[type] += 1
                    break
            #-------------------------
        #-------------------------------
    #----------------------------------
    return countD
#-----------------------------------------

def countLines2(file, typeL):
    countD = {}
    for type in typeL:
        cmd = "zcat "+file+' | grep "'+type+'$" | wc -l'
        countD[type] = int(list(os.popen(cmd))[0])
    #with gzip.open(file, 'rb') as f:
    #    for line in f:
    #        line = line.rstrip()
    #        for type in typeL:
    #            if line.endswith(type):
    #                countD[type] += 1
    #                break
            #-------------------------
        #-------------------------------
    #----------------------------------
    return countD

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    label = options.label
    summary_only = options.summary
    op = options.op
    genome = options.genome
    full = options.full
    per_lb = options.per_lb
    per_ub = options.per_ub
    per_step = options.per_step
    map_reads = options.map_reads
    if map_reads.isdigit():
        map_reads = int(map_reads)
    else:
        for line in open(map_reads):
            if line.startswith("Number of paired-end alignments with a unique best hit"):
                map_reads = int(line.split(':')[1])
                break
    assert isinstance(map_reads, int), map_reads
    verbose = options.verbose
    global debug
    debug = options.debug
    typeL = ['CpG', 'CHH', 'CHG']
    #-----------------------------------
    summary = op + '.summary.xls'
    summary_fh = open(summary, 'w')
    if not summary_only:
        downSample(file, per_lb, per_ub, per_step, op, genome)
    print >>summary_fh, "Per\tvariable\tvalue\tsample"
    print "Per\tvariable\tvalue\tsample"
    for per in range(per_lb, per_ub+per_step, per_step):
        per = per / 100
        prefix = op+'_'+str(per)
        file = prefix+'.SingleCmet.gz'
        print >>sys.stderr, file
        reads = int(per*map_reads)
        countD = countLines(file, typeL)
        for type in typeL:
            print >>summary_fh, "{}\t{}\t{}\t{}".format(reads, type, countD[type], label)
            print "{}\t{}\t{}\t{}".format(reads, type, countD[type], label)

    countD = countLines(full, typeL)
    for type in typeL:
        print >>summary_fh, "{}\t{}\t{}\t{}".format(map_reads, type, countD[type], label)
        print "{}\t{}\t{}\t{}".format(map_reads, type, countD[type], label)
    summary_fh.close()
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


