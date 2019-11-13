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
    This is used to extract overrepresented sequences from the HTML
    output of fastqc.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from bs4 import BeautifulSoup

def fprint(content):
    print json_dumps(content,indent=1)

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file_fastqc.html"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="The HTML output of fastqc.")
    parser.add_option("-f", "--fastqc-data", dest="fastqc_data",
        help="The fastqc_data.txt output by FastQC.")
    parser.add_option("-a", "--adaptor", dest="adaptor",
        default="default", 
        help="A two column file containing the adaptor sequences used \
by fastq. Default: FastQC_dir/Configuration/fastqc_adaptor.fa")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def findComm(adaptor_part, adaptor_full):
    #print "---------------------"
    #print adaptor_part
    #print adaptor_full
    #print "---------------------"
    if adaptor_part.find(adaptor_full) != -1:
        return adaptor_full
    elif adaptor_full.find(adaptor_part) != -1:
        return adaptor_part
    else:
        len_part = len(adaptor_part) - 8
        for i in range(1,len_part):
            if adaptor_full.find(adaptor_part[i:]) != -1:
                return adaptor_part[i:]
    return ''
#--------------findComm---------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    adaptor = options.adaptor
    if adaptor == 'default':
        adaptor = os.path.realpath(os.popen("which fastqc").readlines()[0].strip()).rsplit('/', 1)[0]+"/Configuration/fastqc_adaptor.fa"
        print >>sys.stderr, adaptor
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    adaptorD = dict([line.strip().split('\t') for line in open(adaptor)])
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    fastqc_data = options.fastqc_data
    fh2 = open(fastqc_data)
    line = fh2.readline()
    while 1:
        if line.find(">>Adapter Content") == 0:
            break
        line = fh2.readline()
    adaptorL = set()
    while 1:
        line = fh2.readline()
        lineL = line.strip().split('\t')
        len_lineL = len(lineL)
        if line.find('#Position') == 0:
            headerL = lineL
            continue
        for i in range(1, len_lineL):
            #print >>sys.stderr, lineL[i]
            if float(lineL[i]) > 0.01:
                adaptorL.add(headerL[i])
        if line.find(">>END_MODULE")==0:
            break
    fh2.close()
    #print >>sys.stderr, adaptorL
    for adaptor_name in adaptorL:
        if adaptor_name in adaptorD:
            adaptor = adaptorD[adaptor_name]
            print ">%s\n%s" % (adaptor_name, adaptor)
        else:
            print sys.stderr, "Unexist "+adaptor_name
            sys.exit(1)
    soup = BeautifulSoup(fh, 'html.parser')
    tableL = soup.find_all('table')
    if len(tableL) > 1:
        table = tableL[1]
        aDict = {}
        id = 0
        for tr in table.find_all('tr'):
            id += 1
            aDict[id] = []
            for td in tr.find_all('td'):
                aDict[id].append(td.get_text())
        #-------------END reading file----------
        keyL = aDict.keys()
        keyL.sort()
        for key in keyL:
            valueL = aDict[key]
            if len(valueL) >= 1:
                adaptor_name = valueL[3].split('(')[0].strip()
                if adaptor_name in adaptorD:
                    #adaptor = findComm(valueL[0], adaptorD[adaptor_name])
                    adaptor = adaptorD[adaptor_name]
                    #print adaptor
                    #print '========='
                    if adaptor and len(adaptor) >= 10:
                        print ">%d %s\n%s" % (key, adaptor_name, adaptor)
                        #print '========='
                        #print
                    else:
                        print >>sys.stderr, file, adaptor_name, \
                            valueL[0], adaptorD[adaptor_name]
            #0--------------------
        #--------------------------------
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


