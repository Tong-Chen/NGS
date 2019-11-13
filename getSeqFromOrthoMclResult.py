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
    This is designed to parse orthmcl results.

    Input file format:
    
    cluster_name<colon><any blank>spe1<vertical_line>prot1<any blank>spe2<verticial_line>prot2<any blank>.....

    C10000: Aco|Aco000153.1 Aco|Aco004369.1 Aco|Aco010005.1
    C10001: Aco|Aco000153.1 Cla|Cla004369.1 Dec|Dec010005.1

    Tasks:

    1. Extract protein or nucleotide sequences for given clusters
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
import re
from parseOrthoMclResult import readFastaDir

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
        metavar="FILEIN", help="Output of `orthomclMclToGroups`.")
    parser.add_option("-P", "--directory-prot", dest="dir_prot",
        help="Directory containing all protein sequences used \
for `orthoMcl.sh`. All sequences have a suffix `.fasta`.")
    parser.add_option("-N", "--directory-nucl", dest="dir_nucl",
        help="Directory containing all nucleotide sequences used \
for `orthoMcl.sh`. All sequences have a suffix `.fasta`.")
    parser.add_option("-o", "--output-prefix", dest="outp",
        help="Prefix for output files.")
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
    outp = options.outp
    dir_prot = options.dir_prot
    dir_nucl = options.dir_nucl
    #-----------------------------------------
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    protD = {}
    nuclD = {}
    dataD = []
    if dir_prot:
        outp_fh = open(outp+'.prot.fa', 'w')
        protD = readFastaDir(dir_prot)
        dataD.append([outp_fh, protD])
    if dir_nucl:
        outp_fh = open(outp+'.nucl.fa', 'w')
        nuclD = readFastaDir(dir_nucl)
        dataD.append([outp_fh, nuclD])
    if debug:
        print >>sys.stderr, protD
    #--------------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------

    for line in fh:
        try:
            cluster, protein = line.split(':')
        except ValueError:
            print >>sys.stderr, line
            sys.exit(1)
        proteinL = protein.strip().split()
        for protein in proteinL:
            for fh, seqD in dataD:
                print >>fh, ">%s_%s\n%s" % (cluster, protein, seqD[protein])
    #-------------END reading file----------
    for fh, seqD in dataD:
        fh.close()

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


