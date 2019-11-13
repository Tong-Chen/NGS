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
    This is designed to merge multiple columns file using given map file.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
import re

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
    parser.add_option("-i", "--map-file", dest="filein",
        metavar="FILEIN", help="Map file with first line as header line containing unique labels for each file needed to be merged.")
    parser.add_option("-f", "--mat-file", dest="mat_file",
        help="Multiple columns files (separated by ',' or ' ') with first column as index.")
    parser.add_option("-l", "--file-label", dest="file_label",
        help="Labels for each file (separated by ',' or ' ').")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readInfile(file, label, fileD):
    '''
    fileD = {label: {'header':header_line, 'naline':na_line, key:value}}
    '''
    aDict = {}
    header = 1
    for line in open(file):
        line = line.strip()
        if header:
            headerL = line.split('\t')
            aDict['header'] = '\t'.join(headerL[1:])
            len_h = len(headerL)-1
            labelC = '\t'.join([label]*len_h)
            aDict['label'] = labelC
            naline = ['NA']*len_h
            aDict['naline'] = '\t'.join(naline)
            header -= 1
            continue
        #-------------------------
        key, value = line.strip().split('\t', 1)
        aDict[key] = value
    fileD[label] = aDict
#--------------------------------
def readMapFile(mapfile):
    '''
    mapL = [(label1, id1), (label2, d2), (., .)]
    '''
    mapL = []
    header = 1
    for line in open(mapfile):
        lineL = line.strip().split('\t')
        if header:
            headerL = lineL
            header -= 1
            continue
        tmpMapL = [(key, value) for key, value in zip(headerL, lineL)]
        mapL.append(tmpMapL)
    return headerL, mapL

#--------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    mat_file = options.mat_file
    mat_fileL = re.split(r'[, ]*', mat_file)
    file_label = options.file_label
    file_labelL = re.split(r'[, ]*', file_label)
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    headerL, mapL = readMapFile(file)
    fileD = {}
    '''fileD = {label: {'header':header_line,'label':label_line,'naline':na_line, key:value}}'''
    for mat_file, file_label in zip(mat_fileL, file_labelL):
        readInfile(mat_file, file_label, fileD)
    #-------------------------------------------------    
    print "Label\t{}\tType".format('\t'.join([fileD[i]['label'] for i in headerL]))
    print "ID\t{}\tType".format('\t'.join([fileD[i]['header'] for i in headerL]))
    for tmpL in mapL:
        '''
        mapL = [(label1, id1), (label2, d2), (., .)]
        '''
        tmp_outputL = [tmpL[0][1]]
        labelL = []
        for label, id in tmpL:
            if id not in fileD[label]:
                value = fileD[label]['naline']
            else:
                value = fileD[label][id]
                labelL.append(label)
            #---------------------------------
            tmp_outputL.append(value)
        #----------------------
        tmp_outputL.append(','.join(labelL))
        print '\t'.join(tmp_outputL)
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


