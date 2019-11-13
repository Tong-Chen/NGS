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
    This is specially designed fot ChenWei Data
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from scipy import stats

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
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    ttest = stats.ttest_ind
    header = 1
    for line in fh:
        lineL = line.split()
        len_lineL = len(lineL)
        if header:
            header -= 1
            keyL = lineL[:4]
            for i in range(4, len_lineL, 4):
                typeL = ['FM', 'CD', "Diff", "P_value"]
                sample = lineL[i].rsplit('_', 1)[0]
                keyL.extend([sample+'_'+i for i in typeL])
            print '\t'.join(keyL)
            continue
        #--------------------------------
        keyL = lineL[:4]
        #grp = 
        for i in range(4, len_lineL, 4):
            FM = [float(lineL[i]), float(lineL[i+1])]
            CD = [float(lineL[i+2]), float(lineL[i+3])]
            FM_ave = sum(FM)/2
            CD_ave = sum(CD)/2
            diff = CD_ave - FM_ave
            if diff >= 0.2 or diff <= -0.2:
                pvalue = ttest(FM, CD)[1]
            else:
                pvalue = 1
            value = "%.2f\t%.2f\t%.2f\t%.4f" % (FM_ave, CD_ave, diff, pvalue)
            keyL.append(value)
        #------------------------
        newLine = '\t'.join(keyL)
        print newLine
        #NA_cnt = newLine.count('NA')
        #if NA_cnt + 4 < len_lineL: 
        #    print newLine
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


