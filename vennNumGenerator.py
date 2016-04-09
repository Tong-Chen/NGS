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
    This is designed to generate venn numbers and items for a set of
    samples.

Input file:
    item1   samp1
    item2   samp1
    item3   samp1
    item1   samp2
    item3   samp2
    item4   samp2
    item1   samp3
    item4   samp3
    item5   samp3



'''

from itertools import combinations
import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import re
#from bs4 import BeautifulSoup

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
        metavar="FILEIN", help="A two columns file with \
the first column containing items and the second column \
containing group informations.")
    parser.add_option("-g", "--group", dest="group",
        help="The groups one want to compare. \
<,> or < > separated list. \
<samp1,samp2> will generate samp1_specific, samp2_specific and \
samp1_samp2_common. \
<samp1,samp2,samp3> will generate 7 types of combinations including \
samp1_specific, samp2_specific, samp3_specific, samp1_samp2_common,  \
samp1_samp3_common, samp2_samp3_common, samp1_samp2_samp3_common.")
    parser.add_option("-P", "--plot", dest="plot",
        default=0, help="Plot venn diagram. Default 0 \
means no plotting. Accept 1 to plot. Currently only less \
than 5-sets VennDiagram is supported.")
    parser.add_option("-p", "--output-prefix", dest="prefix",
        help="Specify output prefix")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def difference(groupFirst, groupAll, itemD):
    '''
    Find items common in groupFirst, but not in other groups.
    '''
    groupFirstS = itemD.get(groupFirst[0], set())
    for group in groupFirst[1:]:
        groupFirstS = groupFirstS.intersection(itemD.get(group, set()))

    otherGroupS = set()
    for group in groupAll:
        if group in groupFirst:
            continue
        otherGroupS = otherGroupS.union(itemD.get(group, set()))    
    
    return groupFirstS.difference(otherGroupS)

#----------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    label = ['-a', '-b', '-c', '-d', '-e']
    file = options.filein
    groupL = re.split("[, ]*", options.group.strip())
    plot  = int(options.plot)
    prefix = options.prefix
    output = file + '.' + prefix + '.vennDiagram.xls'
    fh_out = open(output, 'w')
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    itemD = {}
    for line in fh:
        item, group = line.split()
        if group not in itemD:
            itemD[group] = set()
        itemD[group].add(item)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    combineL = []
    len_group = len(groupL)
    for i in range(1, len_group+1):
        combineL.extend(list(combinations(groupL, i)))
    #-------------------------
    for aTuple in combineL:
        itemS = difference(aTuple, groupL[:], itemD)
        print >>fh_out, "%s-specific_to_others\t%d\t%s" % \
            ('_'.join(aTuple), len(itemS), ','.join(itemS))
    fh_out.close()
    #--------------------------------
    if plot and len_group < 6:
        cmd = ['s-plot vennDiagram -f ', file, '-p', prefix]
        for i in range(len_group):
            cmd.append(label[i])
            cmd.append(groupL[i])
        os.system(' '.join(cmd))
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


