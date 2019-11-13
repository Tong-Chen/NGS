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
    This is designed to sort matrix first by the sum of all values and
    then by each row.
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
        metavar="FILEIN", help="A data matrix. \
<-> representing <STDIN>.")
    parser.add_option("-H", "--header", dest="header",
        type="int", default=1, help="Set the number of header \
lines. Default 1.")
    parser.add_option("-F", "--function", dest="func",
        default="sum", help="Given a function with each row as \
input and return value as a sort parameter. Default <sum> to get \
the summary of each row. Other functions are under developing.")
    parser.add_option("-s", "--sort-by-sum", dest="sort_sum",
        default=1, help="Default 1 meaning sort by given function first.\
Accept 0 to sort by others first.")
    parser.add_option("-D", "--digit-type", dest="digit_trans",
        default='int', help="Specify the digit type: \
<int>(default) or <float>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_false", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    header = int(options.header)
    func = options.func
    sort_sum = options.sort_sum
    digit_trans = options.digit_trans
    verbose = options.verbose
    global debug
    debug = options.debug
    funcD = {'sum': sum}
    digit_transD = {'int': int, 'float': float}
    digit_trans = digit_transD[digit_trans]
    func = funcD[func]
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    head_line = ''
    contentL = []
    for line in fh:
        if header:
            head_line = line
            header -= 1
            continue
        lineL = line.split()
        key = lineL[0]
        dataL = [digit_trans(i) for i in lineL[1:]]
        func_v = func(dataL)
        contentL.append([key, dataL, func_v])
    #-------------END reading file----------
    if sort_sum:
        contentL.sort(key=lambda x: (x[2], x[1]))
    else:
        contentL.sort(key=lambda x: (x[1], x[2]))

    #---------Output--------------------------------
    print head_line, 
    for content in contentL:
        dataL = [str(i) for i in content[1]]
        print "%s\t%s" % (content[0], '\t'.join(dataL))

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


