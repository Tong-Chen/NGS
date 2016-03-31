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
    This is designed to extract connections between two group genes by
    giving network.

Input file (At least two columns containing gene pairs are needed.):

RelationNumberOfReferences      Type    Gene1   rela    Gene2
1       Binding SRSF1   ----    ALKBH5
1       Expression      PRMT7   --->    ALKBH5
1       Regulation      ALKBH5  --->    SRPK1
1       Regulation      ALKBH5  --->    SRSF1
1       positive Expression     MAT1A   --+>    METTL3
1       positive Expression     METTL3  --+>    ARNTL
1       positive Expression     METTL3  --+>    PER2
2       Binding Gm14292 ----    METTL3

Output:


'''

import sys
import os
from time import localtime, strftime 
from json import dumps as json_dumps
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

def fprint(content):
    print json_dumps(content, indent=1)

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
        metavar="FILEIN", help="A network file containing at least \
two columns separated by tab with each represeent a list of paired genes.")
    parser.add_option("-U", "--upstream-list", dest="up_file",
        help="The file contains upstream genes .")
    parser.add_option("-D", "--downstream-list", dest="dw_file",
        help="The file contains downstream genes.")
    parser.add_option("-u", "--upstream-col", dest="up_col",
        help="The column contains upstream genes (1-based).")
    parser.add_option("-d", "--downstream-col", dest="dw_col",
        help="The column contains downstream genes (1-based).")
    parser.add_option("-r", "--relationship-col", dest="re_col",
        help="The column contains relationships (1-based). \
If not given, a '--' will be used to represent relationships.")
    parser.add_option("-H", "--header", dest="header",
        default=1, help="The number of header lines to skip, \
default 1.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def traverse(netDict, up_geneD, dw_geneD):
    '''
    traverse(netDict, up_geneD, dw_geneD)
    
    Traverse netDict to screen regulatory relationship between genes
    in up_geneD and dw_geneD.
    '''
    relationD = {}
    for up_gene in up_geneD:
        if up_gene in netDict:
            relationD[up_gene] = []
            up_gene_targetD = netDict[up_gene]
            for dw_gene in up_gene_targetD:
                if dw_gene in dw_geneD:
                    relationD[up_gene].append([up_gene,
                        up_gene_targetD.keys[dw_gene], dw_gene])
                else:

#-------------------------------------------

def readNetwork(fh, header, up_col, dw_col, re_col, relationship='--'):
    netDict = {}
    '''
    1 means Gene1 located upstream Gene2
    -1 means Gene1 located downstream of Gene2
    {
        Gene1: {
            Gene2: [(relation, 1), (relation, -1)]
            Gene3: [(relation, 1)]
        }
        Gene2: {
            Gene1: [(relation, 1), (relation, -1)]
        }
    }
    '''
    for line in fh:
        if header:
            header -= 1
            continue
        #-----------------------------------
        lineL = line.rstrip().split('\t')
        up_gene = lineL[up_col]
        dw_gene = lineL[dw_col]
        if re_col != '':
            relationship = lineL[re_col]
        if up_gene not in netDict:
            netDict[up_gene] = {}
        if dw_gene not in netDict[up_gene]:
            netDict[up_gene][dw_gene] = [(relationship, 1)]
        else:
            tmpL = netDict[up_gene][dw_gene]
            netDict[up_gene][dw_gene] = filter(lambda x:x[1]!=-1, tmpL)
            netDict[up_gene][dw_gene].append((relationship, 1))
        #---------------------------------------------------------
        if dw_gene not in netDict:
            netDict[dw_gene] = {}
        if up_gene not in netDict[dw_gene]:
            netDict[dw_gene][up_gene] = [(relationship, -1)]
        else:
            netDict[dw_gene][up_gene].append((relationship, -1))
        #----------------------------------------------------------
    fprint(netDict)
    return netDict
#------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    up_file = options.up_file
    dw_file = options.dw_file
    up_col = int(options.up_col) - 1
    dw_col = int(options.dw_col) - 1
    re_col = ''
    if options.re_col != None:
        re_col = int(options.re_col) - 1
    header = int(options.header)
    verbose = options.verbose
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    netDict = readNetwork(fh, header, up_col, dw_col, re_col)
    up_geneD = dict([(line.strip(), 1) for line in open(up_file)])
    dw_geneD = dict([(line.strip(), 1) for line in open(dw_file)])
    traverse(netDict, up_geneD, dw_geneD)
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


