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

This is designed to transfer PSL file to bed file for locating query
sequences in genome.

Potential bugs:
    1. Only psl file outputted from Blat is tested. Blat generated psl
    file may have different format as described in UCSC for negative
    strand coordinates. 
    As mentioned in UCSC, "Be aware that the coordinates for a negative
    strand in a PSL line are handled in a special way. In the qStart
    and qEnd fields,  the coordinates indicate the position where the
    query matches from the point of view of the forward strand,  even
    when the match is on the reverse strand. However,  in the qStarts
    list, the coordinates are reversed."
    However, the Blat outputted psl is handled in a normal way. 
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
        metavar="FILEIN", help="The psl file mainly outputed \
from blat.")
    parser.add_option("-s", "--skip-lines", dest="header",
        default=5, help="The number of lines you want to skip \
before real processing. Default 5 since blat output with \
headers in 5 lines.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def transferPSLtoBed(filehandle, header):
    '''
    See here for PSL format: 
    http://genome.ucsc.edu/FAQ/FAQformat.html#format2
    '''
    aDict = {}
    for line in filehandle:
        if header:
            header -= 1
            continue
        #---------------------
        #using '\t' in case there are gaps in query name
        lineL = line.strip().split('\t')
        strand = lineL[8]
        name = lineL[9]
        qSize = lineL[10]
        chr = lineL[13]
        blockSize = [int(i) for i in lineL[18].strip(',').split(',')]
        tStart = [int(i) for i in lineL[20].strip(',').split(',')]
        count = len(tStart)
        for i in range(count):
            tmpL = '\t'.join([chr, str(tStart[i]),
                str(tStart[i]+blockSize[i]), name,
                qSize, strand])
            if name in aDict:
                aDict[name].append(tmpL)
            else:
                aDict[name] = [tmpL]
        #------------------------------------------
    return aDict
#-----------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    header = int(options.header)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = transferPSLtoBed(fh, header)
    for key, valueL in aDict.items():
        print '\n'.join(valueL)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
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



