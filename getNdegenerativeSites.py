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
    This is designed to get 4-fold degenerate sites (4DTv) or codons
    for given CDS sequences.
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
        metavar="FILEIN", help="A multiple aignment FASTA file")
    parser.add_option("-m", "--multiple-aligned-fasta", dest="msa",
        default='phylip',  help="Specify the type of given \
multiple alignment file.Currently <phylip>(default) or <fasta> \
is acccepted. ")
    parser.add_option("-g", "--gap-keep", dest="gap_keep",
        default=1, type='int',  help="Keep sites with gaps (default).\
Accept <0> to delete sites with gaps.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

#def geneDegenerateCodons():
    
#-------------------------------------
def readPhylip(file):
    '''
    Read MSA in phylip format
    '''
    seqD = {}
    fh = open(file)
    contextL = fh.readlines()
    header = contextL[0]
    #print header
    count, length = header.split()
    count = int(count)
    step  = count + 1

    len_context = len(contextL)
    
    for j in range(0, step):
        if j:
            name, seq = contextL[j].strip().split(' ', 1)
            assert name not in seqD, "Duplicate %s" % name
            seqD[name] = ''
            seqL = [seq]
        for i in range(j+step, len_context, step):
            line = contextL[i].strip()
            if i % step == 0:
                assert not line, line 
                continue
            name_ind = i % step
            seqL.append(line)
            #----------------------------
        #----------------------------
        if j:
            seqD[name] = ''.join(seqL).replace(' ', '')
    #----------------------------------
    return seqD
#---------readPhylip-----------------

def readFasta(file):
    '''
    Read MSA in Fasta format
    '''
    seqD = {}
    for line in open(file):
        if line[0] == '>':
            key = line[1:-1]
            assert key not in seqD, "Duplicate seqname %s" % key
            seqD[key] = []
        else:
            seqD[key].append(line.strip())
    #----------------------------------------
    for name, seqL in seqD.items():
        seqD[name] = ''.join(seqL).replace(' ', '')
    return seqD
#------------readFasta--------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    msa  = options.msa
    msaD = {'phylip':readPhylip, 'fasta':readFasta}
    gap_keep = options.gap_keep
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    degenerateCodons = {'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 
                        'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 
                        'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 
                        'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 
                        'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 
                        'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 
                        'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 
                        'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 
                        '---':'#'
                        }
    #-----------------------------------
    
    seqD = msaD[msa](file)
    #print seqD

    nameL  = seqD.keys()
    nameL.sort()
    valueL = seqD.values()

    len_seq = len(valueL[0])
    
    savePos = []

    for i in range(0,len_seq,3):
        print >>sys.stderr, "Position %d" % i
        codonS = set()
        for seq in valueL:
            codon = seq[i:i+3]
            aa    = degenerateCodons.get(codon, '*') 
            codonS.add(aa)
            if debug:
                print >>sys.stderr, "\tcodon:%s aa:%s" % (codon, aa)
        #-----------------------------
        codonL = list(codonS)
        if '*' in codonL:
            continue
        if '#' in codonL:
            if gap_keep:
                codonL.remove('#')
            else:
                continue
        if len(codonL) == 1:
            savePos.append(i+2)
            print >>sys.stderr, "\n\tSaved\n"
    #------------------------------------------
    for name in nameL:
        print ">%s" % name
        seq = seqD[name]
        tmpL = [seq[pos] for pos in savePos]
        print ''.join(tmpL)

    if debug:
        print >>sys.stderr, savePos
        start = 0
        tmpL = []
        for pos in savePos:
            for i in range(start, pos):
                tmpL.append('-')
            tmpL.append('N')
            start = pos+1
        for i in range(start, len_seq):
            tmpL.append('-')
        print >>sys.stderr, '>sumamry'
        print >>sys.stderr, ''.join(tmpL)
        for name in nameL:
            seq = seqD[name]
            print >>sys.stderr, '>%s' % name
            print >>sys.stderr, seq
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


