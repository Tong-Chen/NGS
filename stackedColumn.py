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
Functionla description

This is designed to generate data for plotting stacked column.

***Past test 2013-11-22***
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

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
        metavar="FILEIN", help="The file want to compute. \
Multiple files separated by <,> are accepted.")
    parser.add_option("-f", "--format", dest="format",
        metavar="FASTA", default='FATSA', help="The format of files \
given. Default is one-line-FASTA. For FATSA file, \
all sequences belong to one item must in one \
line. Currently only FATSA and plain file is allowed. \
Actually FASTA is treated as plain file without checking the \
legal format.")
    parser.add_option("-t", "--type", dest="type",
        metavar="COUNT", default='count', help="The type of data to be outputed. \
<count> represents raw count of each element. \
<percentage> represents the relative count of each element.")
    parser.add_option("-F", "--output-format", dest="of",
        default='matrix', help="Default output normal matrix-table. \
if <melt> given, a melted matrix will be outputed.")
    parser.add_option("-s", "--Set", dest="Set",
        default=0, help="Only used when -F is melt. A word to \
represent the attribute of your data. If multiple files are \
supplied to -i, related string for <Set> is needed.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def stackColumnsFromFasta(file):
    aDict = {}
    for line in open(file):
        if line[0] == '>':
            continue
        else:
            count = 0
            for i in line[:-1]:
                count += 1
                if count in aDict:
                    if i in aDict[count]:
                        aDict[count][i] +=1
                    else:
                        aDict[count][i] = 1
                else:
                    aDict[count] = {}
                    aDict[count][i] = 1
        #----------------------------------
    #-------END reading----------
    return aDict
#-------------------------------------------

def outputStackColumns(aDict, type, of, Set, noHead):
    posKeyL = aDict.keys()
    posKeyL.sort() 
    elementL = [j for i in posKeyL for j in
        aDict[i].keys()]
    elementL = list(set(elementL))
    elementL.sort()
    if of == 'matrix':
        if noHead == 0:
            print "pos\t%s" % '\t'.join(elementL)
        for pos in posKeyL:
            tmpL = [aDict[pos].get(ele,0) for ele in elementL]
            if type == 'count':
                tmp = '\t'.join([str(count) for count in tmpL])
            elif type == 'percentage':
                sumL = sum(tmpL) * 1.0
                tmp = '\t'.join([str(i/sumL) for i in tmpL])
            #-----------------------------------------
            print "%s\t%s" % (str(pos), tmp)
        #---------------------------------------
    elif of == 'melt':
        lenEle = len(elementL)
        if noHead == 0:
            print "pos\tvariable\tvalue\tSet" 
        for pos in posKeyL:
            tmpL = [aDict[pos].get(ele,0) for ele in elementL]
            if type == 'count':
                tmpL = [str(i) for i in tmpL]
            elif type == 'percentage':
                sumL = sum(tmpL) * 1.0
                tmpL = [str(i/sumL) for i in tmpL]
            #-----------------------------------------
            for i in range(lenEle):
                print "%s\t%s\t%s\t%s" % \
                    (str(pos),elementL[i], tmpL[i], Set)
    #--------------------------------------------
#------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    fileL = options.filein.split(',')
    format = options.format
    type = options.type
    of = options.of
    if options.Set != None:
        Set = options.Set.split(',')
        assert len(Set) == len(fileL)
    else:
        Set = ''
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    #--------------------------------
    noHead = 0
    for file in fileL:
        aDict = stackColumnsFromFasta(file)
        if Set != '':
            outputStackColumns(aDict, type, of, Set[noHead], noHead)
        else:
            outputStackColumns(aDict, type, of, '', noHead)
        #--------------------------------
        noHead += 1
    #-------------END reading file----------
    #----close file handle for files-----
    #-----------end close fh-----------
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



