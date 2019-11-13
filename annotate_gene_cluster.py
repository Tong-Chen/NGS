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
    This is designed to annotate gene cluster file generated by find_gene_cluster.sh.

Cluster file

>Cluster 10
MRNA_000218_3   MRNA_000949_1   MRNA_001063_1   MRNA_000953_1
>Cluster 18
MRNA_000399_1   MRNA_000415_2   MRNA_000214_1   MRNA_000422_1
>Cluster 28
MRNA_000538_1   MRNA_000744_2   MRNA_000742_1
>Cluster 31
MRNA_000371_1   MRNA_000137_2   MRNA_000466_2   MRNA_000226_1
>Cluster 38
MRNA_000437_1   MRNA_000284_1   MRNA_000574_1
>Cluster 39
MRNA_000667_1   MRNA_000549_3   MRNA_000517_1   MRNA_000799_2
>Cluster 40
MRNA_000591_3   MRNA_000442_1   MRNA_001167_5   MRNA_000460_1

Annotation file may contain multiple columns separated by <tab>.
The first column matches with gene names (MRNA_000218_3) in each cluster, 
other columns contains different levels of annotation.
Clusters will be omitted if no gene is annotated.

ID  Annotate    Pfam
MRNA_000218_3   CYP *****
MRNA_000218_4   CYP *****
MRNA_000218_5   CYP *****
MRNA_000218_6   CYP *****
MRNA_000218_7   CYP *****

'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0


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
        metavar="FILEIN", help="Cluster file with format showed above.")
    parser.add_option("-a", "--annotation-file", dest="annotation_file",
        help="Annotation file with format specified above")
    parser.add_option("-n", "--column-num", dest="column",
        default=0, 
        help="Specify which column would be used for annotation. Default all columns. Accept a series of numbers separated by ',' like <2,5,7> meaning 2nd, 3rd and fifth column would be used for annotation. Default all columns would be used.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

#-----------------------------------
def readAnno(annotation_file, column, header=1):
    annoD = {}
    
    keepL = []
    for line in open(annotation_file):
        lineL = line.strip().split('\t')
        if header:
            #headerD = dict(enumerate(lineL))
            headerL = lineL[:]
            header -= 1
            len_col = len(headerL)
            if column != 0:
                columnL = [int(i)-1 for i in column.split(',')]
            else:
                columnL = [i for i in range(1, len_col)]
            keepL = [headerL[i] for i in columnL]
            continue
        #----------------------------------------------------
        key = lineL[0]
        assert key not in annoD, "Duplicate "+ key
        annoD[key] = {}
        for i in columnL:
            secKey = headerL[i]
            secValue = lineL[i]
            if secValue in ['.', '-'] or len(secValue)==1:
                secValue = 'NA'
            annoD[key][secKey] = secValue
    #-------------------------------------
    return annoD, keepL
#-----------------------------------



def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    annotation_file = options.annotation_file
    column = options.column
    
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    """
    annoD = {'gene1':{'anno_typ1':'anno_value1', 'anno_typ2':'anno_value2'}}
    keepL = ['anno_typ1', 'anno_typ2']
    """

    annoD, keepL = readAnno(annotation_file, column)
    #-----------end close fh-----------

    for line in open(file):
        if line[0] == '>':
            key = line[1:-1]
        else:
            geneL = line.split()
            len_geneL = len(geneL)
            annoL = [annoD.get(gene, 'NA') for gene in geneL]
            if len_geneL > 1 and annoL.count('NA') < len_geneL:
                print ">"+key
                print line,
                for key in keepL:
                    tmpL = []
                    for subD in annoL:
                        if subD == 'NA':
                            tmpL.append('NA')
                        else:
                            tmpL.append(subD[key])
                    print key+':'+'\t'.join(tmpL)
                    #print '\t'.join([suD[key] for subD in annoL])
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


