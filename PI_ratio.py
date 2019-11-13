#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright 2018, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    this is used to compute PI ratio output by vcftools.

CHROM   BIN_START   BIN_END N_VARIANTS  PI
E_H535995_c3_g4 1   50000   12  8.37972e-05
E_H535995_c3_g3 1   50000   5   1.75147e-05
E_H535995_c3_g1 1   50000   49  0.00020649
E_H544449_c0_g6 1   50000   2   7.40741e-07

Group_file
#divider    dividen type
Cultivation_type    Wilf vcftools

1. Header line is optional (lines start with '#' will be skipped)
2. Column order is fixed as header line
3. if type==vcftools, i=SNP_analysis/Cultivation_type.vcftools.population_statistics.windowed.pi, j=SNP_analysis/Wild.vcftools.population_statistics.windowed.pi, m=Cultivation_type, n=Wild

'''

import sys
import os
import re
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import pandas as pd
import numpy as np
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
    parser.add_option("-f", "--group-file", dest="group_file",
        metavar="FILEIN", help="A 3 column file as described above. This will overwrite i, j, m, n parameters.")
    parser.add_option("-i", "--divider-file", dest="divider_file",
        metavar="FILEIN", help="A list of files separated by blank or comma used as divider.")
    parser.add_option("-j", "--dividen-file", dest="dividen_file",
        metavar="FILEIN", help="A list of files separated by blank or comma used as dividen.")
    parser.add_option("-m", "--divider-name", dest="divider_name",
        metavar="FILEIN", help="A list of sample names separated by blank or comma used as divider names in file. Please use legally simple words.")
    parser.add_option("-n", "--dividen-name", dest="dividen_name",
        metavar="FILEIN", help="A list of files separated by blank or comma used as dividen.")
    parser.add_option("-x", "--divide-all", dest="divide_all",
        default=False, action="store_true", help="Default one divider divides on dividen correspondingly. If set <divide-all>,  one divide will divide all dividen one by one.")
    parser.add_option("-o", "--output-prefix", dest="prefix",
        default="PI_ratio", help="Output prefix")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    if not options.group_file:
        assert options.divider_file != None, "A filename needed for -i"
        assert options.dividen_file != None, "A filename needed for -j"
        assert options.divider_name != None, "A filename needed for -m"
        assert options.dividen_name != None, "A filename needed for -n"
    return (options, args)
#--------------------------------------------------------------------

def PI_ratio(divider_file,  divider_name, dividen_file, dividen_name, prefix):
    divider_matrix = pd.read_table(divider_file, header=0, index_col=[0, 1, 2], usecols=[0, 1, 2, 4])
    divider_matrix.columns = [divider_name]
    dividen_matrix = pd.read_table(dividen_file, header=0, index_col=[0, 1, 2], usecols=[0, 1, 2, 4])
    dividen_matrix.columns = [dividen_name]

    matrix = divider_matrix.join([dividen_matrix],  how="outer")
    matrix = matrix.fillna(0)
    matrix[divider_name+'.divide.'+dividen_name] = matrix[divider_name]/matrix[dividen_name]
    matrix = matrix[~matrix.isin([np.nan, np.inf, -np.inf]).any(1)]
    output = '.'.join([prefix, divider_name, 'divide', dividen_name, "txt"])
    matrix.to_csv(output,  sep=b"\t")

#--------------------------------------------------------------------
def parseGroupFile(group_file):
    divider_fileL = []
    dividen_fileL = []
    divider_nameL = []
    dividen_nameL = []
    for line in open(group_file):
        if line[0] == '#':
            continue
        divider, dividen, type1 = line.strip().split()
        if type1 == 'vcftools':
            divider_file = ['SNP_analysis/', divider, '.vcftools.population_statistics.windowed.pi']
            divider_fileL.append(''.join(divider_file))
            divider_nameL.append(divider)
            dividen_file = ['SNP_analysis/', dividen, '.vcftools.population_statistics.windowed.pi']
            dividen_fileL.append(''.join(dividen_file))
            dividen_nameL.append(dividen)
        else:
            print >>sys.stderr, "Unsupported type "+type1
        return divider_fileL, dividen_fileL, divider_nameL, dividen_nameL 
#--------------------------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    group_file = options.group_file
    if group_file:
        divider_fileL, dividen_fileL, divider_nameL, dividen_nameL = parseGroupFile(group_file)
    else:
        divider_fileL = re.split(r'[, ]*',options.divider_file)
        dividen_fileL = re.split(r'[, ]*',options.dividen_file)
        divider_nameL = re.split(r'[, ]*',options.divider_name)
        dividen_nameL = re.split(r'[, ]*',options.dividen_name)
    #print divider_fileL, dividen_fileL, divider_nameL, dividen_nameL 
    divide_all    = options.divide_all
    prefix        = options.prefix
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    file_len = len(divider_fileL)
    if divide_all:
        for i in range(file_len):
            divider_file = divider_fileL[i]
            divider_name = divider_nameL[i]
            for j in range(file_len):
                dividen_file = dividen_fileL[j]
                dividen_name = dividen_nameL[j]
                PI_ratio(divider_file,  divider_name, dividen_file, dividen_name, prefix)
    else:
        for i in range(file_len):
            divider_file = divider_fileL[i]
            divider_name = divider_nameL[i]
            dividen_file = dividen_fileL[i]
            dividen_name = dividen_nameL[i]
            PI_ratio(divider_file,  divider_name, dividen_file, dividen_name, prefix)
    #-----------end close fh-----------
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
    if verbose:
        #print("--Successful %s" % strftime(timeformat, localtime()), file=sys.stderr)
        print >>sys.stderr, "--Successful %s" % strftime(timeformat, localtime())

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    with open('python.log', 'a') as fh:
        #print ("%s\n\tRun time : %s - %s " % \
        #(' '.join(sys.argv), startTime, endTime), file=fh)
        print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


