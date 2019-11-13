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
    Designed for one time usage.
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
        metavar="FILEIN", help="")
    parser.add_option("-s", "--sum-only", dest="sum_only",
        default=False, action="store_true", metavar="SUM-ONLY", 
        help="Specify to onply do summary.")
    parser.add_option("-c", "--column-anno", dest="col_anno",
        help="Column anno file")
    #parser.add_option("-n", "--number", dest="number",
    #    type="int", help="Supply an int number")
    #parser.add_option("-c", "--choice", dest="choice",
    #    type="choice", choices=["a", "b", "c"], 
    #    default="a", help="Supply an int number")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    sum_only = options.sum_only
    anno = options.col_anno
    op   = file
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
    dataD = {}
    attribD = {}
    """ 
    dataD = {symbol:{samp1: {pos1:1, pos2:2, }, 
                     samp2:{pos1:1, pos2:2}, 
                     samp3:{pos1:1, pos2:2}, }
            }
    attribD = {symbol:{'strand':'+', 'pos':[pos1, pos2, pos3, .]}}
    """
    len_lineL = 0 
    for line in fh:
        lineL = line.strip().split('\t')
        if not len_lineL:
            len_lineL = len(lineL)
        if header:
            headerL = lineL
            headerD = dict(enumerate(headerL))
            header -= 1
            continue
        #------------------------------------
        sym, strand, chr, pos, _ = lineL[0].split('@', 4)
        pos = int(pos)
        if sym not in dataD:
            dataD[sym] = {}
            for samp in headerL[1:]:
                dataD[sym][samp] = {}
            attribD[sym] = {'strand': strand, 'pos':[pos]}
        else:
            attribD[sym]['pos'].append(pos)
        for i in range(1, len_lineL):
            samp = headerD[i]
            dataD[sym][samp][pos] = lineL[i]
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    sumD = {}
    
    for symbol, attribSubD in attribD.items():
        #sumD[symbol] = {}
        if not sum_only:
            out_each_file = op +'.'+symbol+'.xls'
            out_each = open(out_each_file, 'w')
        strand = attribSubD['strand']
        posL   = attribSubD['pos']
        if strand == '-':
            posL.sort(reverse=True)
        else:
            posL.sort()
        if not sum_only:
            print >>out_each, "%s\t%s" % ("Samp", '\t'.join([str(i) for i in posL]))
        for samp in headerL[1:]:
            if samp not in sumD:
                sumD[samp] = {}
            sumD[samp][symbol] = sum([float(i) for i in dataD[symbol][samp].values()])
            if not sum_only:
                tmpL = [dataD[symbol][samp][pos] for pos in posL]
                print >>out_each, "%s\t%s" % (samp, '\t'.join(tmpL))
        if not sum_only:
            out_each.close()
            cmd = ['s-plot pheatmap -A 45 -R FALSE -E png -f', out_each_file, '-t', symbol]
            os.system(' '.join(cmd))
    ###--------multi-process------------------
    out_file = op+'.totalSum.xls'
    out = open(out_file, 'w')
    symbolL = attribD.keys()
    print >>out, "samp\t%s" % '\t'.join(symbolL)
    '''
    sumD = {'samp1':[syn1:1, syn2:2, ]}
    '''
    for samp in headerL[1:]:
        tmpL = [str(sumD[samp][i]) for i in symbolL]
        print >>out, "%s\t%s" % (samp, '\t'.join(tmpL))
    out.close()
    cmd = ["transpose.py", out_file, ">"+out_file+'.tmp']
    os.system(' '.join(cmd))
    cmd = ["/bin/mv", out_file+'.tmp', out_file]
    os.system(' '.join(cmd))

    cmd = ['s-plot pheatmap -a TRUE -b FALSE -R TRUE -H FALSE -E png -d row -u 15 -v 15 -A 45 -f', out_file, '-Q', anno]
    os.system(' '.join(cmd))
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


