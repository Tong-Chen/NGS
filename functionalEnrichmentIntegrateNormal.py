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
    This is designed to integrate the functional enrichment result of
    functionalEnrichmentAnalysis.
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
    parser.add_option("-a", "--all_de_file", dest="allDE",
        default="allDE", help="The all.DE file generated \
using DESeq.sh.")
    parser.add_option("-g", "--gotxt", dest="gotxt",
        default="/MPATHB/soft/homer/data/GO/GO.txt", 
        help="A list of annotation files supplied by homer.\
Default is GO.txt file in <data> directory of <homer>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.allDE != None, "A filename needed for -a"
    return (options, args)
#--------------------------------------------------------------------

def readGOtxt(gotxt, header=1):
    GOL = []
    for line in open(gotxt):
        if header:
            header -= 1
            continue
        lineL = line.split('\t')
        #file, descp, link, linkN = lineL[2:6]
        GOL.append(lineL[2:6][:])
    #--------------------
    return GOL
#-------------------------

def mergeResult(anno, comparePL, go_dir, allDE, top=50, max_len=70):
    fileL = ['/'.join([go_dir, i, anno]) for i in comparePL]
    i = -1
    unanno = -1
    annoL = ['Term\tneg_log10FDR\tCount\tSample']
    for file in fileL:
        i += 1
        if os.path.isfile(file):
            header = 1
            count = 0
            for line in open(file):
                if header:
                    header -= 1
                    continue
                #----------------------------
                lineL = line.split('\t')
                tmpLine = '\t'.join([lineL[1][:max_len], lineL[3], 
                    lineL[5], comparePL[i]])
                annoL.append(tmpLine)
                count += 1
                if count > top: break
            #----------END reading one file-----------
        else:
            annoL.append('\t'.join(["No enrichment", '0', '0',
                comparePL[i]]))
            unanno += 1
    #---------------------------------------------------------
    if unanno < i and len(annoL) > 2:
        file96 = allDE+"." +anno
        fh = open(allDE+"."+anno, 'w')
        print >>fh, '\n'.join(annoL)
        fh.close()
        height = len(annoL) / 3
        if height < 8:
            height = 10
        if height < 25:
            height = 25
        height = str(height)
        width = str(30)
        if len(fileL) == 1:
            dv_plot = ["s-plot scatterplotDoubleVariable -f", file96, 
                "-o Sample -v Term -c neg_log10FDR -s Count -w 25 -a", height, 
                "-E pdf"]
        else:
            dv_plot = ["s-plot scatterplotDoubleVariable -f", file96, 
                "-o Sample -v Term -c neg_log10FDR -s Count -w 30 -a", height, 
                "-E pdf -R 90 -H 0 -V 1"]
        print ' '.join(dv_plot)
        os.system(' '.join(dv_plot))
        convert = ["convert -density 150 -quality 90",
            file96+".scatterplot.dv.pdf", file96+".scatterplot.dv.png"]
        os.system(' '.join(convert))
#-------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    allDE = options.allDE
    de_dir = os.path.dirname(allDE)
    gotxt = options.gotxt
    if gotxt:
        GOL = readGOtxt(gotxt) 
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------

    comparePL = []
    for line in open(allDE):
        compareP = line.strip().split()[1]
        if compareP not in comparePL:
            comparePL.append(compareP)
   #-------------------------------------

    go_dir = de_dir + "/go/"
    for funcL in GOL:
        mergeResult(funcL[0], comparePL, go_dir, allDE)

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


