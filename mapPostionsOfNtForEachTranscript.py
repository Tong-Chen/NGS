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

    This program is designed to generate the coordinates for each
    nucleotide of given mRNAs.
FASTA format:
    >NM_001195662 gene=Rp1 CDS=55-909
    AAGCTCAGCCTTTGCTCAGATTCTCCTCTTGATGAAACAAAGGGATTTCTGCACATGCTTGAGAAATTGC
    AGGTCTCACCCAAAATGAGTGACACACCTTCTACTAGTTTCTCCATGATTCATCTGACTTCTGAAGGTCA
    AGTTCCTTCCCCTCGCCATTCAAATATCACTCATCCTGTAGTGGCTAAACGCATCAGTTTCTATAAGAGT
    (strings between '>' and first space ' ' will be used as sequence
    name; This name should match the names in Bed12 file. 
    Numbers after 'CDS=' would be used as additional indicator
    if sequence name is duplicated in bed12 file.)

Bed12 format:
    Standard format
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
        metavar="FASTA", help="Transcript sequence in FASTA format \
with their name as real names. Usually generated using <gffread> from \
<cufflinks>.")
    parser.add_option("-g", "--gtf-file", dest="gtf",
        metavar="GTF", help="The gtf file used to generate the fasta \
sequence. Only 'exon' information is used.")
    parser.add_option("-b", "--bed12", dest="bed12",
        metavar="BED12", help="The bed12 file. If you use gtf file \
to generate the fasta sequence, please use < \
gtfToGenePred mm9.gene.refgene.gtf test.pred;  \
genePredToBed test.pred test.bed> to get Bed12 from GTF.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readFasta(fh):
    seqDict = {}
    for line in fh:
        if line[0] == '>':
            locus = line.strip()[1:]
        else:
            if locus not in seqDict:
                seqDict[locus] = [line.strip()]
            else:
                seqDict[locus].append(line.strip())
        #------------------------------------------
    #--------------------------------------------
    for key, valueL in seqDict.items():
        seqDict[key] = ''.join(valueL)
    return seqDict
#--------------------------------------------

def readGTF(gtf):
    '''
    unfinished
    '''
    for line in open(gtf):
        lineL = line.split('\t')
        chr = lineL[0]
        start = lineL[3]
        end = lineL[4]
        strand = lineL[6]
        attributeL = lineL[8].split('"')
        gene = attributeL[1]
        transcript = attributeL[3]
        locus = transcript + ' gene=' + gene
#----------------------------------------------

def readBed12(bed12):
    '''
    Return a dict in format like
    aDict = {tr:
                {cds1: 
                    [   chr,
                        strand, 
                        ntCnt,
                        [(exon_start,exon_end), 
                         (          ,        )
                        ]
                    ], 
                cds2:[]
                }
            }
    '''
    aDict = {}
    for line in open(bed12):
        lineL = line.strip().split('\t')
        chr = lineL[0]
        key = lineL[3]
        strand = lineL[5]
        tcStart = int(lineL[1])
        tcEnd   = int(lineL[2])
        tsStart = int(lineL[6])
        tsEnd   = int(lineL[7])
        exonCnt = int(lineL[9])
        exonSize = [int(i) for i in lineL[10].strip(',').split(',')]
        ntCount = sum(exonSize)
        exonStart = [int(i) for i in lineL[11].strip(',').split(',')]
        exonL = [(tcStart+exonStart[i], tcStart+exonStart[i]+exonSize[i])
                for i in range(exonCnt)]
        if strand == '+':
            cds_startpos = str(tsStart - tcStart + 1)
            cds_endpos   = str(tsEnd - tcEnd + 1) 
        elif strand == '-':
            cds_startpos = str(tcEnd - tsEnd + 1)
            cds_endpos   = str(tcStart - tsStart + 1) 
        #------------------------------------------
        secondKey = cds_startpos + '-' + cds_endpos
        if key not in aDict:
            aDict[key] = {}
        if secondKey not in aDict[key]:
            aDict[key][secondKey] = [chr,strand, ntCount, exonL]
        else:
            print >>sys.stderr, "There may be something wrong, \
since duplicated key %s found with transcript %s" % (secondKey, key)
            sys.exit(1)
    #-----------------------------------
    return aDict
#-------------------------------------

def mapSeqToCoord(seqDict, coordDict):
    '''
    '''
    posDict = {}
    for name,seq in seqDict.items():
        firstSpace = name.find(' ')
        if firstSpace != -1:
            tr = name[:firstSpace]
            cds_p = name.find('CDS=')
            if cds_p != -1:
               cds = name[cds_p+4:]
            else:
               cds = ''
        else:
            tr = name
            cds = ''
        #---------------------------
        innD = coordDict[tr]
        if len(innD) == 1:
            posDict[name] = mapSeqToCoordEach(name,len(seq),innD.values()[0])
        else:
            if cds in innD:
                posDict[name] = mapSeqToCoordEach(name,len(seq),innD[cds])
            else:
                print "Ambiguous sequence", tr
        #--------END for--------------
    return posDict
#-------------------------------------------

def mapSeqToCoordEach(name,len_seq,coordL):
    chr     = coordL[0]
    strand  = coordL[1]
    ntCount = coordL[2]
    exonposL  = coordL[3]
    assert len_seq == ntCount, name
    posL = [chr,strand]
    tmpL = []
    for i in exonposL:
        tmpL.extend([str(j) for j in range(i[0], i[1])])
    if strand == '-':
        tmpL.reverse()
    posL.extend(tmpL)
    return posL
#------------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    gtf  = options.gtf
    bed12 = options.bed12
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    seqDict = readFasta(fh)
    if bed12:
        coordDict = readBed12(bed12)
    #--------------------------------------------
    posDict = mapSeqToCoord(seqDict, coordDict)
    for name,coord in posDict.items():
        print '>%s' % name
        print ','.join(coord)
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



