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
    This is designed to classify P450 gene to family and clan.
    
    The family is defined by the alphabets (CYP) and numbers before
    the first alphabet. 

    The clans are defined according to one plant journal, 
    and summarized in file
    "/MPATHB/resource/P450/Plant_CYPclan_CYPfamily.txt".

Family  Gene
CYP71   CYP71A1
CYP749  CYP749A22
CYP714  CYP714C2
CYP82   CYP82A3
CYP78   CYP78A3
CYP71   CYP71A2

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from bs4 import BeautifulSoup
import re

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
        metavar="FILEIN", help="The file with one column \
containing CYP gene information.")
    parser.add_option("-c", "--column", dest="col",
        help="The column containing CYP genes. 1-based.")
    parser.add_option("-s", "--column-str", dest="column_str",
        help="A string to specify the column containing CYP genes. \
Normally <Gene_name> for my own usage.")
    parser.add_option("-C", "--clan", dest="clan",
        default="/MPATHB/resource/P450/Plant_CYPclan_CYPfamily.txt", 
        help="CYP450 clan and family file.")
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
    if options.col:
        col  = int(options.col) - 1
    elif options.column_str:
        col_str = options.column_str
        col = -1
    else:
        assert 1==0, "-c or -s should be specified"

    clan = options.clan
    pat  = re.compile(r"CYP[^A-Z]*", flags=re.IGNORECASE)
    clanD = {}
    header = 1
    for line in open(clan):
        if header:
            header -= 1
            continue
        #----------------
        clan, familyL = line.split()
        familyL = familyL.split(';')
        for family in familyL:
            if family not in clanD:
                clanD[family] = clan
            else:
                print >>sys.stderr, "Ambigious clan for %s" % family
                sys.exit(1)
    #------------------------------
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    header = 1
    for line in fh:
        if col != -1:
            lineL = line.strip().split('\t', col+1)
        if header:
            if col == -1:
                lineL = line.split('\t')
                count = 0
                for ele in lineL:
                    if ele == col_str:
                        col = count
                        count += 1
                        break
                assert count, "Unrecgonizable {} in {}".format(col_str, line) 
            #------------------------------------
            lineL.insert(col, "CYP_family")               
            lineL.insert(col, "CYP_clan")               
            print '\t'.join(lineL)
            header -= 1
            continue
        #---------------------------------
        cyp = lineL[col].upper()
        lineL[col] = cyp
        family = pat.match(cyp)
        if family:
            family = family.group()
        else:
            family = cyp
        clan = clanD.get(family, family)
        lineL.insert(col, family)               
        lineL.insert(col, clan)               
        print '\t'.join(lineL)
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


