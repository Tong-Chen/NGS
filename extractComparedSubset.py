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
    This is designed to extract compared groups information from all groups for substream analysis.

    The output wolud be a json file.


data_matrix
frag_len        variable        value   set
100     B59_20G 1119    Fragment length distrib (count)
101     B59_20G 2439    Fragment length distrib (count)
102     B59_20G 3487    Fragment length distrib (count)
103     B59_30G 3394    Fragment length distrib (count)
104     B59_30G 3991    Fragment length distrib (count)
105     B59_30G 3624    Fragment length distrib (count)
106     F59_20G 3779    Fragment length distrib (count)
107     F59_20G 3847    Fragment length distrib (count)
108     F59_20G 4115    Fragment length distrib (count)
109     F59_20G 4286    Fragment length distrib (count)
110     F59_20G 4590    Fragment length distrib (count)

sampleFile (normally two columns is OK. If there are nested group info, extra columns can be added.)

Samp    Grp1    Grp2
B59_20G   B59   20G
B59_30G   B59   30G
F59_20G   F59   20G
F59_30G   F59   30G

compare_pair
B59 F59
20G 30G


norm_factor
A   1000
B   1000
C   1000
D   1000
E   1000
F   1000
G   1000
H   1000

Output JSON

[
    {"grp1": "grp1_matrix"}, 
    {"grp2": "grp2_matrix"}, 
    {"grp3": "grp3_matrix"}
]

'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from math import log
from statsmodels.stats.multitest import multipletests
from scipy import stats

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
    usages = "\n%prog -i file\nzcat gzfile | %prog -i -"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A data matrix. Gzipped file should be zcat to STDIN.")
    parser.add_option("-k", "--input-key-column", dest="input_key",
        help="Key column to match with sampleFile. 1-based.")
    parser.add_option("-S", "--input-split-column", dest="input_split",
        help="Split column to separate into multiple files. 1-based")
    parser.add_option("-s", "--sample-file", dest="sampleFile",
        metavar="Samp-cond", help="Format as described above to specify group information for each sample. ")
    parser.add_option("-c", "--compare-file", dest="compareFile",
        metavar="compare-pair", help="Format as described above to specify which two samples needed to be compared. If no compareFile specified, each two groups will be compared.")
    parser.add_option("-o", "--output-prefix", dest="op",
        help="Output prefix for all files.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the programe")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def readSampleFile(sampleFile, sampleL=[]):
    '''
    Samp    GrpA    GrpB
    A   Grp1    grp1
    B   Grp2    grp1
    C   Grp1    grp1
    D   Grp1    grp2
    E   Grp2    grp2
    F   Grp3    grp3
    G   Grp4    grp3
    H   Grp4    grp4
    '''
    sampleD = {}
    if sampleFile:
        header = 1
        for line in open(sampleFile):
            if header:
                header -= 1
                continue
            lineL = line.split()
            samp = lineL[0]
            grpL = lineL[1:]
            for grp in grpL:
                if grp not in sampleD:
                    sampleD[grp] = [samp]
                else:
                    sampleD[grp].append(samp)
            #-------------------
        #-------END for-------------
    else:
        sampleD = dict([i, [i]] for i in sampleL)
    return sampleD
#---------------------------------------------
def readCom_pair(compareFile, grpL=[]):
    '''
    Grp1    Grp2
    Grp3    Grp4
    '''
    if compareFile:
        return [tuple(line.split()) for line in open(compareFile)]
    else:
        grpL.sort()
        len_grpL = len(grpL)
        compareL = []
        for i in range(len_grpL-1):
            for j in range(i+1, len_grpL):
                compareL.append((grpL[i], grpL[j]))
    return compareL
#---------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file         = options.filein
    input_key    = int(options.input_key)-1
    input_split    = options.input_split
    if input_split:
        input_split = int(input_split)-1
    else:
        input_split = -1
    sampleFile   = options.sampleFile
    compareFile  = options.compareFile
    verbose      = options.verbose
    op           = options.op
    global debug
    debug = options.debug
    #-----------------------------------
    #--------------------------------
    samp2posD = {}
    '''
    sampleD = {'grp1':[A, B], 'grp2':[C, D]}
    compareL = [(grp1, grp2), (grp1, grp3)]
    normD = {'A':1, 'B':1, 'C':1}
    compareD = {(grp1, grp2):{id1:[[1, 2], [3, 4]], id2:[[1, 3], [5, 6]]}}
    '''
    sampleD  = readSampleFile(sampleFile)
    compareL = readCom_pair(compareFile)
    if debug:
        print >>sys.stderr, "sampleD", sampleD
        print >>sys.stderr, "compareL", compareL
    header = 1
    fileD = {}
    '''
    fileD = {'ID1': [line1, line2, line3, ...], 
             'ID2': [line1, line2, line3, ...]
            }
    '''
    for line in open(file):
        line  = line.strip()
        if header:
            header -= 1
            headerLine = line
            continue
        lineL = line.split('\t')
        ID    = lineL[input_key]
        if ID not in fileD:
            fileD[ID] = []
        fileD[ID].append(line)
    #-----------------------------------
    json_sum = open(op+".json", 'w')
    opL = []
    for grpL in compareL:
        opD = {}
        grp_name = '.'.join(grpL)
        tmp_file = op + '._.' + grp_name + '.xls'
        opD[grp_name] = tmp_file
        tmp_file_fh = open(tmp_file, 'w')
        print >>tmp_file_fh, headerLine+"\tGroup_extra"
        for grp in grpL:
            for samp in sampleD[grp]:
                for line in fileD[samp]:
                    print >>tmp_file_fh, line+"\t"+grp
        tmp_file_fh.close()
        opL.append(opD)
    #------------------------------
    json_dump(opL, json_sum, indent=2)
    json_sum.close()
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


