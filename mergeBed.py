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
'''
Functionla description
Merge bed regions together for binned regions of mRNA.
Those regions spearated by introns will be identified and added
together.

Binned-regions(At least 4 columns, more than 6 columns are allowed.
The last column will be taken as a number which will be averaged in
merging process.)

chr1	84034117	84034142	NM_001003948_21304.UTR3__69__51@897504  0	-
chr1    84034142        84034167    NM_001003948_21304.UTR3__69__52@897505  0       -

Template-bed
chr5    31351012        31351045        NM_027855_3.UTR5        0   +
chr5    31356740        31356996        NM_027855_3.UTR3        0   +
chr5    130695613       130695789       NM_001081394_4.UTR5     0   +
chr5    130705318       130705337       NM_001081394_4.UTR5     0   +
chr5    130717168       130719635       NM_001081394_4.UTR3     0   +

'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    desc = ""
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="Binned-regions", help="Binned regions. \
'-' as STDIN accepted.")
    parser.add_option("-t", "--template-file", dest="template",
        metavar="Template-bed", help="template file")
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
    template = options.template
    verbose = options.verbose
    debug = options.debug
    #--------get templateD--------------
    temD = {}
    for line in open(template):
        lineL = line.split()
        nameL = lineL[3].split('.')
        gene  = nameL[0]
        if gene not in temD:
            temD[gene] = []
        temD[gene].append((int(lineL[1]), int(lineL[2])))
    #--------generate junctionD------------
    if verbose:
        print >>sys.stderr, "temD"
        print >>sys.stderr, temD
    juncD = {}
    for gene, valueL in temD.items():
        valueL.sort(key=lambda x: x[0])
        if verbose:
            print >>sys.stderr, 'valueL'
            print >>sys.stderr, valueL
        juncD[gene] = []
        start = valueL[0][1]
        end   = ""
        for posS in valueL[1:]:
            end = posS[0]
            #assert start < end, "%d %d %s " % (start, end, gene)
            if start < end:
                juncD[gene].append((start, end))
                start = posS[1]
            if start == end:
                start = posS[1]
    #---------get binned regions------------
    if verbose:
        print >>sys.stderr, "juncD"
        print >>sys.stderr, juncD
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #---generate binD------------------------
    binD = {}
    for line in fh:
        lineL = line.split()
        name  = lineL[3].split('.',1)
        gene  = name[0]
        key = (int(lineL[1]), int(lineL[2]))
        if gene not in binD:
            binD[gene] = {}
        #--------------------------------------
        binD[gene][key] = lineL
    #-------------END reading file----------
    #------------begin merge---------------
    for gene, valueD in binD.items():
        mergeD = {}
        juncDL = juncD[gene]
        valueD_k = valueD.keys()
        valueD_k.sort(key=lambda x: x[0])
        initial_k = list(valueD_k[0])
        mergeD[tuple(initial_k)] = [valueD[tuple(initial_k)]]
        for followK in valueD_k[1:]:
            if initial_k[1] == followK[0]:
                tmpInitial_k = tuple(initial_k[:])
                initial_k[1] = followK[1]
                mergeD[tuple(initial_k)] = [i116 for i116 in mergeD[tmpInitial_k]]
                mergeD[tuple(initial_k)].append(valueD[followK])
                mergeD.pop(tmpInitial_k)
            else:
                potentialJunc = (initial_k[1], followK[0])
                if potentialJunc in juncDL:
                    tmpInitial_k = tuple(initial_k[:])
                    initial_k[1] = followK[1]
                    mergeD[tuple(initial_k)] = [i116 for i116 in mergeD[tmpInitial_k]]
                    mergeD[tuple(initial_k)].append(valueD[followK])
                    mergeD.pop(tmpInitial_k)
                else:
                    initial_k = list(followK)
                    mergeD[tuple(initial_k)] = [valueD[followK]]
            #--------------------------------------------------
        #----END merge one gene-------------------------------------
        if verbose:
            print >>sys.stderr, "mergeD"
            print >>sys.stderr, mergeD
        #-------------output final merged results---------------
        i = 0
        for keyS, valueL in mergeD.items():
            i += 1
            j = 0
            initialL = valueL[0]
            lenInitialL = len(initialL)
            if lenInitialL > 6:
                lenInitialL = 6
            nameL = initialL[3].split('.', 1)
            name = nameL[0]
            type = set([nameL[1].split('__')[0]])
            mean = [float(initialL[-1])]
            peakL = []
            fullLen = 0
            for followL in valueL[1:]:
                if initialL[2] == followL[1]:
                    initialL[2] = followL[2]
                    type.add(followL[3].split('.', 1)[1].split('__')[0])
                    mean.append(float(followL[-1]))
                else:
                    j += 1
                    initialL[3] = "__".join([name,str(i),str(j)]) 
                    initialL[4] = str(sum(mean)*1.0/len(mean))
                    #print '\t'.join(initialL[:6])
                    peakL.append(initialL[:lenInitialL])
                    fullLen += int(initialL[2]) - int(initialL[1])
                    initialL = followL
                    #nameL = initialL[3].split('.', 1)
                    #name = nameL[0]
                    #type = [nameL[1].split('__')[0]]
                    type.add(followL[3].split('.', 1)[1].split('__')[0])
                    mean = [float(initialL[-1])]
            #---------------------------------------------------------------
            j += 1
            initialL[3] = "__".join([name,str(i),str(j)]) 
            initialL[4] = str(sum(mean)*1.0/len(mean))
            #print '\t'.join(initialL[:6])
            peakL.append(initialL[:lenInitialL])
            fullLen += int(initialL[2]) - int(initialL[1])
            type = list(type)
            type.sort()
            type = '-'.join(type)
            for splitPeak in peakL:
                print "%s\t%s\t%d" % ('\t'.join(splitPeak),
                    type, fullLen) 
         #------------------------END output----------------------- 
    #---END merge--------------------------------------

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



