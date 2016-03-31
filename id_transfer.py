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
    This is first designed to transfer GTF id to entrez id.

It needs two input file:
1. ID file 
#--No header lines needed-------#
#--The first column is the ID column, will be used only--
#--Other columns will be output as input-----------------
#--File content------------------
AT1G04770.1     cbf0._vs_.cbf3_up
AT1G05660.1     cbf0._vs_.cbf3_up
AT1G07350.2     cbf0._vs_.cbf3_up
AT1G10150.1     cbf0._vs_.cbf3_up
AT1G10640.1     cbf0._vs_.cbf3_up
AT1G11700.1     cbf0._vs_.cbf3_up
AT1G12805.1     cbf0._vs_.cbf3_up
AT1G14200.1     cbf0._vs_.cbf3_up
AT1G14430.1     cbf0._vs_.cbf3_up
AT1G15825.1     cbf0._vs_.cbf3_up
AT1G15830.1     cbf0._vs_.cbf3_up
AT1G19670.1     cbf0._vs_.cbf3_up
*************OR******************
#--File content------------------
AT1G04770.1
AT1G05660.1
AT1G07350.2
AT1G10150.1
AT1G10640.1
AT1G11700.1
AT1G12805.1
AT1G14200.1
AT1G14430.1
AT1G15825.1
AT1G15830.1
AT1G19670.1
AT1G22640.1

2. ID map file
#--One Header line needed------------------#
#--At least two columns are required---#
#--The first column should match ID column in ID file--#
#--The second column will be used for the transfer   --#
ID   EntrezGene
AT3G18710   821402
AT4G25880   828694
AT4G25880   828694
AT4G25880   828694
AT1G71695   843498

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
    parser.add_option("-i", "--id-file", dest="filein",
        metavar="ID", help="The ID file with format specified above.")
    parser.add_option("-m", "--id-map", dest="idmap",
        metavar="ID-MAP", help="The ID-map file with format specified above.")
    parser.add_option("-c", "--id-col", dest="id_col",
        default=2, help="Specify the column of target IDs. Default <2> \
indicating the second column.")
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
    idmap = options.idmap
    id_col = int(options.id_col) - 1
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    unmatchL = set()

    #idDict = dict([[line.strip().split('\t')[0], \
    #                line.strip().split('\t')[id_col]] \
    #                for line in open(idmap)])
    idDict = {}
    for line in open(idmap):
        lineL = line.strip('\n').split('\t')
        key = lineL[0]
        map = lineL[id_col]
        if map:
            if key in idDict:
                idDict[key].add(map)
            else:
                idDict[key] = set([map])
    #----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    for line in fh:
        lineL = line.strip().split("\t", 1)
        id = lineL[0]
        idNew = idDict.get(id, "Unmatch")
        #print id, idNew
        if idNew == "Unmatch" or idNew == "":
            unmatchL.add(line)
        else:
            for id in idNew:
                lineL[0] = id
                print '\t'.join(lineL)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    unmatchL = list(unmatchL)
    if unmatchL:
        unmatch = open(file + '.unmatch', 'w')
        print >>unmatch, ''.join(unmatchL)
        unmatch.close()
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


