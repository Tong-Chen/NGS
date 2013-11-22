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

This is designed to transfer a mulicolumn file into a matrix using
given string as head line. A matching will be given 1 and other 0.

Input file:
NM_001001144_24108      NM_001001144    Scap    Coding_exon
NM_001001182_4562       NM_001001182    Baz2b   Coding_exon
NM_001001183_27504      NM_001001183    Tmem204 Coding_exon
NM_001001295_8082       NM_001001295    Dis3l   Coding_exon
NM_001001297_27162      NM_001001297    Thnsl1  Coding_exon
NM_001001321_25094      NM_001001321    Slc35d2 TsTS
NM_001001327_28966      NM_001001327    Vkorc1l1        TsTS
NM_001001333_13296      NM_001001333    Hexdc   TsTS

Output file
Gene    TcSS    UTR5    Coding_exon TsTS    UTR3    
NM_001001144_24108  0   0   1   0   0
NM_001001182_4562   0   0   1   0   0

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
        metavar="FILEIN", help="Format as description, \
multiple-column file without header.")
    parser.add_option("-s", "--string", dest="head",
        default="TcSS;UTR5;Coding_exon;TsTS;UTR3", 
        metavar="str1;str2;str3", 
        help="A semicolon separated file in a order-mattered way. Default \
'TcSS;UTR5;Coding_exon;TsTS;UTR3'.")
    parser.add_option("-m", "--match-column", dest="mc",
        default=4, help="The column used for matching \
the string given to -s. 1-based, default 4.")
    parser.add_option("-n", "--name-column", dest="nc",
        default=1, help="The column for output as the name column. \
Default 1 (1-based).")
    parser.add_option("-F", "--file-contain-full-name-list", dest="ff",
        help="A file containing full list of names specific by \
name-column. This is unnecessary. However,  if this is given, \
it will be used as a template for the name column,  <umv> would be given \
to all names not exist in --file. < - > represents STDIN is accepted.")
    parser.add_option("-a", "--value-for-matching", dest="mv",
        default='1', help="The value to represent matching. \
Default 1 (string format).")
    parser.add_option("-b", "--value-for-unmatch", dest="umv",
        default='0', help="The value to represent matching. \
Default 0 (string format).")
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
    ff = options.ff
    verbose = options.verbose
    debug = options.debug
    headL = options.head.split(';')
    mc = int(options.mc) - 1
    nc = int(options.nc) - 1
    mv = options.mv
    umv = options.umv
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {}
    for line in fh:
        lineL = line.strip().split('\t')
        name = lineL[nc]
        match = lineL[mc]
        if name not in aDict:
            aDict[name] = [match]
        else:
            aDict[name].append(match)
        #-----------------------------------
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    #------------Begin output the ones not in --file-----------------
    print 'head\t%s' % '\t'.join(headL)
    if ff:
        if ff == '-':
            fh = sys.stdin
        else:
            fh = open(ff)
        tmpV = '\t'.join([umv for i in headL])
        for line in fh:
            key = line.strip()
            if key not in aDict:
                print '%s\t%s' % (key, tmpV)
            else:
                tmpL = [key]
                valueL = aDict[key]
                for i in headL:
                    if i in valueL:
                        tmpL.append(mv)
                    else:
                        tmpL.append(umv)
                #-------------END for--------
                print '\t'.join(tmpL)
            #--------------------------------
        #----------------------------------------
        if ff != '-':
            fh.close()
    #------------Begin output the ones in --file-----------------
    else:
        keyL = aDict.keys()
        keyL.sort()
        for key in keyL:
            tmpV = [key]
            valueL = aDict[key]
            for i in headL:
                if i in valueL:
                    tmpV.append(mv)
                else:
                    tmpV.append(umv)
            #-------------END for--------
            print '\t'.join(tmpV)
        #----------------------------------
    #-----------------------------------------------------
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



