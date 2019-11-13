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
    This is designed to transfer an expression matrix for boxplot
    with set information added.

    A matrix to 3 or 4 coumns file.

---SampleFile-------------------
Samp    conditions
T0_1    T0
T0_2    T0
T0_3    T0
T2_1    T2
T2_2    T2
T2_3    T2

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
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
        metavar="FILEIN", help="An expression matrix file.")
    parser.add_option("-s", "--sampleFile", dest="sample",
        metavar="sampleFile", 
        help="sampleFile with format specified above. Optional. Only first two columns would be used.")
    parser.add_option("-g", "--geneFile", dest="gene",
        metavar="gene", 
        help="geneFile with one row containing one gene. Optional.")
    parser.add_option("-G", "--geneString", dest="geneS",
        help="Blank (< >) separated a list of genes for extraction. Optional. ")
    parser.add_option("-a", "--addIDCol", dest="addIDCol",
        default=False, action="store_true", 
        help="Add a series of number to the first column (used for importing to mysql).")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_false", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    if options.sample:
        sampD = dict([line.split()[:2] for line in open(options.sample)])
    else:
        sampD = {}
    gene = options.gene
    addIDCol = options.addIDCol
    geneD = {}
    if gene:
        geneD = dict([(line.strip(), 1) for line in open(gene)])
    geneS = options.geneS
    if geneS:
        geneD = dict([(i, 1) for i in geneS.split()])
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    header = 1

    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    header_outputL = []
    id = 0
    if addIDCol:
        header_outputL.append('ID')
    header_outputL.append('gene')
    header_outputL.append('sample')
    header_outputL.append('value')

    if sampD:
        header_outputL.append('Set')

    print '\t'.join(header_outputL)
    for line in fh:
        if header:
            headerL = line.rstrip().split('\t')
            #print headerL
            if debug:
                print >>sys.stderr, headerL
            len_head = len(headerL)
            if debug:
                print >>sys.stderr, len_head
            header -= 1
            continue
        #--------------------------
        lineL = line.split()
        #print lineL
        gene = lineL[0]
        if geneD and gene not in geneD:
            continue
        for i in range(1, len_head):
            samp = headerL[i]
            id += 1
            tmpL = []
            if addIDCol:
                tmpL.append(str(id))
            if debug:
                print >>sys.stderr, lineL
            tmpL.extend([gene, samp, lineL[i]])
            if sampD:
                set = sampD[samp]
                tmpL.append(set)
            print "\t".join(tmpL)
    #-------------END reading file----------
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


