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
    This is designed to check the matching status for given position
    of given miRNAs. 

    It depends on the predicted secondary structure of nafold (.det
    file), a parsed miRNA.dat file to map miRNAs to its
    pre-miRNA(generayed by parsemirBase.dat.embl.py).

Input file:
    See the end of script.

Score formula:
    # match = 2 mismatch = 0 GU = 1

***Test passed***
'''

import sys
import re
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
    parser.add_option("-n", "--nafold-det", dest="det",
        metavar="NAFOLD.DET", help="The det file outputed from nafold. \
Normaly you can run nafold through wrapped nafold.pl.")
    parser.add_option("-m", "--map-file", dest="map",
        metavar="DAT.MAP", help="The parsed miRNA.dat file though \
parsemirBase.py. This is used to map miRNAs to pre-miRNAs.") 
    parser.add_option("-s", "--subset-file", dest="sub",
        metavar="SUB", help="A file containing only part of \
miRNAs.[Unsupported]")
    parser.add_option("-a", "--all-pos-status", dest="all_status",
        metavar="TRUE", help="Return the pairing status of all \
positions. Default FALSE")
    parser.add_option("-p", "--position", dest="pos",
        metavar="POS", help="The positions you want to check. \
Accept discrete numbers like '1,2,4,6' or a range '1-5' or both. \
All mentioned numbers and internal numbers will be included and \
treated as 1-based from 5' end of miRNAs.") 
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.det != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def complement(base1, base2):
    aDict = {'A':'U', 'U':'A', 'C':'G', 'G':'C'}
    if base1 == aDict[base2]:
        return 2
    elif (base1 == 'G' and base2 == 'U') or \
            (base1 == 'U' and base2 == 'G'):
            return 1
    else:
        return 0
#------------end complement---------------
def debug(aDict):
    for name, itemD in aDict.items():
        print name
        keyL = itemD.keys()
        keyL.sort()
        len_keyL = len(keyL) 
        half = len_keyL/2
        assert len_keyL % 2 == 0
        for i in range(half):
            key1 = keyL[i]
            key2 = keyL[len_keyL-i-1]
            print "%d\t--%d--\t%d\t--%d-" % \
                (key1, itemD[key1], key2, itemD[key2])
#------------END debug----------


def readDet(det, pos=1):
    aDict = {}
    pat = re.compile(r'losing pair is ([A-Z])\( *([0-9]*)\)-([A-Z])\( *([0-9]*)\)')
    fh = open(det)
    line = fh.readline()
    while 1:
        while not line.startswith('Structure    1'):
            line = fh.readline()
            if not line:
                break
        if not line:
            break
        line = fh.readline() # ignore the blank line
        #Get the name of the structure
        name = fh.readline().split(' ')[pos] 
        aDict[name] = {} # initialize a dict of a dict
        line = fh.readline() # ignore dG line
        line = fh.readline() # ignore blank line
        #--begin parsing structure--------------
        line = fh.readline() # the first structure line
        while line.strip():  # end structure when meet blank line
            match_obt = pat.search(line)
            if not (line.startswith('Helix') or \
                    line.startswith('External loop') or \
                    line.startswith(' ')):
                assert match_obt, line
            if match_obt:
                base1 = match_obt.group(1)
                pos1  = int(match_obt.group(2))
                base2 = match_obt.group(3)
                pos2  = int(match_obt.group(4))
                # match = 2 mismatch = 0 GU = 1
                assert pos1 not in aDict[name]
                assert pos2 not in aDict[name]
                aDict[name][pos1] = complement(base1, base2)
                aDict[name][pos2] = aDict[name][pos1]
            else: # only for testing
                assert line.find('(') == -1, line
                #print >>sys.stderr, line,
            #-----------------------------
            line = fh.readline()
            if not line:
                break
        #----------END one structure--------------------------
        #debug(aDict)
        if not line:
            break
    #----------------END all structures-------------
    fh.close()
    return aDict
#------------readDet--------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    det = options.det
    map = options.map
    pos = options.pos
    #parse pos
    posL = []
    for i in pos.split(','):
        if i.find('-') != -1:
            start, end = i.split('-')
            start = int(start)
            end   = int(end) + 1
            for i_71 in range(start, end):
                posL.append(i_71)
        else:
            posL.append(int(i))
        #----------------------------
    posL.sort()
    #------------------------------------
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    aDict = readDet(det)
    #--------------------
    print "miR\t%s" % '\t'.join([str(i) for i in posL])
    for line in open(map):
        un,ac,mir,un,pos = line.split()
        start, end = pos.split('..')
        start = int(start)
        #centralL = [i+start for i in posL]
        valueL = [str(aDict[ac].get(start+i-1, 0)) for i in posL]
        print "%s\t%s" % (mir, '\t'.join(valueL))
    #-----------------------------------
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


'''
<hairpin.fa.det>

Structure    1

mmu-mir-504 MI0005515 Mus musculus miR-504 stem-lo
 dG =    -44.30

External loop:      ddG =  -1.40   5 ss bases &  1 closing helices.
Stack:           ddG =  -2.10 External closing pair is C(     3)-G(    76)
Stack:           ddG =  -2.10 External closing pair is U(     4)-A(    75)
Stack:           ddG =  -2.20 External closing pair is G(     5)-C(    74)
Helix:           ddG =  -6.40   4 base pairs.
Bulge loop:      ddG =   2.90 External closing pair is U(     6)-A(    73)
Stack:           ddG =  -2.10 External closing pair is U(     7)-A(    71)
Stack:           ddG =  -3.30 External closing pair is G(     8)-C(    70)
Stack:           ddG =  -3.30 External closing pair is G(     9)-C(    69)
Stack:           ddG =  -2.40 External closing pair is G(    10)-C(    68)
Stack:           ddG =  -0.60 External closing pair is A(    11)-U(    67)
Stack:           ddG =  -1.30 External closing pair is G(    12)-U(    66)
Stack:           ddG =  -2.20 External closing pair is A(    13)-U(    65)
Stack:           ddG =  -3.30 External closing pair is C(    14)-G(    64)
Stack:           ddG =  -3.30 External closing pair is C(    15)-G(    63)
Stack:           ddG =  -2.10 External closing pair is C(    16)-G(    62)
Stack:           ddG =  -2.10 External closing pair is U(    17)-A(    61)
Helix:           ddG = -26.00  12 base pairs.
Interior loop:   ddG =  -1.00 External closing pair is G(    18)-C(    60)
Stack:           ddG =  -1.50 External closing pair is U(    20)-G(    58)
Stack:           ddG =  -2.10 External closing pair is C(    21)-G(    57)
Stack:           ddG =  -2.10 External closing pair is U(    22)-A(    56)
Stack:           ddG =  -3.40 External closing pair is G(    23)-C(    55)
Helix:           ddG =  -9.10   5 base pairs.
Interior loop:   ddG =   0.40 External closing pair is C(    24)-G(    54)
Stack:           ddG =  -2.10 External closing pair is C(    26)-G(    52)
Stack:           ddG =  -2.40 External closing pair is U(    27)-A(    51)
Stack:           ddG =  -2.10 External closing pair is C(    28)-G(    50)
Helix:           ddG =  -6.60   4 base pairs.
Interior loop:   ddG =   3.10 External closing pair is U(    29)-G(    49)
Stack:           ddG =  -2.40 External closing pair is U(    31)-A(    46)
Helix:           ddG =  -2.40   2 base pairs.
Interior loop:   ddG =   0.40 External closing pair is C(    32)-G(    45)
Stack:           ddG =  -2.20 External closing pair is G(    34)-C(    43)
Stack:           ddG =  -1.30 External closing pair is U(    35)-A(    42)
Helix:           ddG =  -3.50   3 base pairs.
Hairpin loop:       ddG =   5.30          Closing pair is A(    36)-U(    41)


<dat.table,  only column 2,3,5 is used>

mmu-let-7g	MI0000137	mmu-let-7g-5p	MIMAT0000121	7..28
mmu-let-7g	MI0000137	mmu-let-7g-3p	MIMAT0004519	63..84
mmu-let-7i	MI0000138	mmu-let-7i-5p	MIMAT0000122	6..27
mmu-let-7i	MI0000138	mmu-let-7i-3p	MIMAT0004520	62..83
mmu-mir-1a-1	MI0000139	mmu-miR-1a-1-5p	MIMAT0016979	10..32
mmu-mir-1a-1	MI0000139	mmu-miR-1a-3p	MIMAT0000123	49..70
mmu-mir-15b	MI0000140	mmu-miR-15b-5p	MIMAT0000124	4..25
mmu-mir-15b	MI0000140	mmu-miR-15b-3p	MIMAT0004521	42..63
mmu-mir-23b	MI0000141	mmu-miR-23b-5p	MIMAT0016980	9..29
mmu-mir-23b	MI0000141	mmu-miR-23b-3p	MIMAT0000125	46..66
'''
