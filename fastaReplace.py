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
    this is designed to replace special amino acids by positions.

Input file:

    Normal FASTA

    >pESCp4_1
    MDSFSLLAALFFISAATWF

    or 

    aligned FASTA

    >pESCp4_1
    MDSFSLLAALFFISAATWF---ISSRRRR----NLPPGPFPYPIVGNMLQLG

Replace file:

1. Currently this type of file is supported.

   * The first line contains two strings (one corresponds to FASTA sequence name, the other correponds to data label) separated by a blank or a tab. This line will be used as name of replaced FASTA sequence.
   * The second line will be ignored.
   * Only the first, third and forth column of other lines will be used. 
     * First column represents amino acid position in original Fasta sequence (1-based, '-' ignored). Only the first appearance of the position will be kept for duplicated positions.
     * The third column represents the amino acid in given position of original Fasta sequence. Only used to make sure the correctness of the program.
     * The forth column represents the new name of correpsoning amino acid. This will appear in new sequences. 

---------------------------
pESCp4_1 ferru
Res_position;Res_name;Res_name;distance
111;MET;M;1
479;PHE;F;1
111;MET;M;2
302;THR;T;2
306;GLU;E;2
362;PRO;P;2
366;LEU;L;2
---------------------------

2. 

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import re
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
        metavar="FILEIN", help="One or multiple sequences Fasta file.")
    parser.add_option("-f", "--replace-file", dest="replace_file",
        help="One or several files separated by <,> or < >.")
    parser.add_option("-t", "--replace-file-type", dest="replace_file_type",
        default=1, help="Currently only one type for replace files is supported. Default 1.")
    parser.add_option("-k", "--keep-original", dest="keep_original",
        default=False, action="store_true", 
        help="Default replace; If tunned on, only As matched to given position will be kept.")
    parser.add_option("-o", "--output-file", dest="output_file",
        help="Name for the output file.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readInReplace_file(replace_fileL, replace_file_type, nameD, posD):
    '''
    P4 ferru
    Res_position;Res_name;Res_name;distance
    111;MET;M;1
    479;PHE;F;1
    111;MET;M;2
    302;THR;T;2
    '''
    for replace_file in replace_fileL:
        if debug:
            print >>sys.stderr, "Reading in "+replace_file
        header = 2
        for line in open(replace_file):
            if header == 2:
                fasta_name, label = line.strip().split()
                new_name = '_'.join([fasta_name, label])
                posD[new_name] = {}
                if fasta_name not in nameD:
                    nameD[fasta_name] = [new_name]
                else:
                    nameD[fasta_name].append(new_name)
                header -= 1
                continue
            elif header == 1:
                header -= 1
                continue
            line = line.strip()
            if line:
                lineL = line.split(';')
                pos = int(lineL[0]) - 1
                if pos not in posD[new_name]:
                    posD[new_name][pos] = lineL[2:]
        #----------------------------
    #--------------------------------------
#--------------------------------------------------------

def replaceSeq(seq, replaceD, keep_original=False):
    '''
    seq = 'VKRRADVYFGRLLALIEGYLNDRIQSRKANPDAPKKDDFLETLVDILNSN'
    replaceD = {'0':['M', '1'], '1':['D', '1'], }
    '''
    pos = 0
    newSeqL = []
    for i in seq:
        if i == '-':
            newSeqL.append(i)
            continue
        replaceL = replaceD.get(pos, '')
        if debug:
            print >>sys.stderr, pos
            print >>sys.stderr, i
            print >>sys.stderr, replaceL
        if replaceL:
            len_replaceL = len(replaceL)
            if len_replaceL == 2:
                assert i == replaceL[0], 'Unmatched amino acid pos={} original={} new={}'.format(pos, i, replaceL[0])
            if keep_original:
                newSeqL.append(i)
            else:
                newSeqL.append(replaceL[-1])
        else:
            newSeqL.append('-')
        pos += 1
    return ''.join(newSeqL)
#''''''''''''''''''''''''''''''''''''''''''''''



def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file          = options.filein
    replace_file  = options.replace_file.strip()
    replace_fileL = re.split(r'[, ]*', replace_file)
    replace_file_type = options.replace_file_type
    keep_original     = options.keep_original
    output_file = options.output_file
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    nameD = {}
    posD = {}
    if debug:
        print >>sys.stderr, replace_fileL
    readInReplace_file(replace_fileL, replace_file_type, nameD, posD)
    if debug:
        print >>sys.stderr, nameD
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    seqL = []
    output_fh = open(output_file, 'w')
    for line in fh:
        if line[0] == '>':
            if seqL:
                seq = ''.join(seqL)
                print >>output_fh, '\n'.join(['>'+key, seq])
                replace_keyL = nameD.get(key)
                if debug:
                    print >>sys.stderr, key, replace_keyL
                if replace_keyL:
                    for replace_key in replace_keyL:
                        replaceD = posD.get(replace_key)
                        if replaceD:
                            newSeq = replaceSeq(seq, replaceD, keep_original)
                            print >>output_fh, '\n'.join(['>'+replace_key, newSeq])
            key = line.strip()[1:]
            seqL = []
        else:
            seqL.append(line.strip())
    #-------------END reading file----------
    if seqL:
        seq = ''.join(seqL)
        print >>output_fh, '\n'.join(['>'+key, seq])
        replace_keyL = nameD.get(key)
        if debug:
            print >>sys.stderr, key, replace_keyL
        if replace_keyL:
            for replace_key in replace_keyL:
                replaceD = posD.get(replace_key)
                if replaceD:
                    newSeq = replaceSeq(seq, replaceD, keep_original)
                    print >>output_fh, '\n'.join(['>'+replace_key, newSeq])
            #-----------------------------------------
        #-----------------------------------

    #----close file handle for files-----
    if file != '-':
        fh.close()
    output_fh.close()
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


