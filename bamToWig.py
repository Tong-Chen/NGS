#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================

import sys
import re
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

def cmdparameter(argv, options, args):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    desc = "This is written to transfer BAM file to \
wig file. It can deal with single-end bam and pair-end bam and \
strand-specific bam and the combination of them. For pair-end bam, a \
GTF file is needed for RNA-Seq. The program will try to assign reads for middle \
regions."
    usages = "%prog [-i SAM file] [-o output/stdout]"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A SAM file with or without header,  or \
- means STDIN. File must be sorted by chromosome")
    parser.add_option("-o", "--output-file", dest="fileout",
        metavar="FILEOUT", help="If not given, STDOUT is used")
    parser.add_option("-g", "--gtf", dest="gtf",
        metavar="GTF", help="When -t is PE and -e is positive, this \
GTF should be supplied. It canbe standard GTF downlaoded from UCSC. \
However, the one outputted by you RNA-Seq would be better (Only with \
expressed transcripts is preferred). File must be sorted by \
chromosome.")
    parser.add_option("-n", "--nucleotide-type", dest="nt",
        metavar="RNA/DNA", default="RNA", 
        help="DNA means ChIP-Seq, RNA means RNA-Seq")
    parser.add_option("-t", "--seq-type", dest="seq_Type",
        metavar="PAIREND/SINGLEEND", default="PE", 
        help="PE for pair-end reads and SE for single-end reads.")
    parser.add_option("-s", "--strand", dest="strand",
        default=True, action="store_false", 
        help="True for strand-specific, False for strandless sequenceing")
    parser.add_option("-e", "--extend", dest="extend", default=0,
        ,type='int', help="A positive number means extending reads to given length \
for SE data or filling in the blank between two PE reads assisted by \
GTF. 0 means no extending or filling. For PE reads,  any positive \
number can be used to represent extend.")
    parser.add_option("-c", "--chrom-size", dest="chromSize",
        help="If -t is SE and -e larger than 0 and no header in SAM \
file. This should be given. One can use \
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
'select chrom, size from mm10.chromInfo' > mm10.genome \
to extract chromosome size.")
    (options, args) = parser.parse_args(argv)
#--------------------------------------------------------------------

def computeRegion(start, cigarL, name): 
    # cigarL:  [('20','M'),('1000','N'),('80','M')]
    regionL = []
    len = 0
    for i in cigarL:
        if i[1] == 'N': #begin settle counts
            assert len > 0, "Unexpected cigarL %s " % name
            regionL.append([start, start+len])
            start = start + len + int(i[0])
            len = 0
        else:
            if i[1] == 'M' or i[1] == 'D':
                len += int(i[0])
            elif i[1] == 'I':
                pass # no adding length
            else:
                print >>sys.stdout, "Unconsidered cigars %s" % name
        #------------------------------------
    #--------------END for-------------------
    #------If no 'N' or deal with the part after the last 'N'--
    regionL.append([start, start+len])
    return regionL
    #--------------------------------------------------------------------
#--------END computeRegion-------------------------------------------------

def computeWigDict(wigDict, pairL):
    #wigDict = {pos:{+:[+,+_e], '-':[ -,-_e]}}
    #pairL   = [[chr,[[start,end],...],xs], [chr,[[start,end],...],xs]]
    
#----------END computeWigDict---------------
def extendWigDict(wigDict,gtf_fh):
    pass

#--------NED extendWigDict---------------------

def main():
    options = {}
    args = []
    cmdparameter(sys.argv[1:], options, args)
    #-----------------------------------
    cigarP = re.compile('([0-9]+)([A-Z])')
    file = options.filein
    output = options.fileout
    readsType = opions.seq_Type
    strand = options.strand
    extend = options.extend
    gtf = options.gtf
    nt = options.nt
    cs = options.chromSize
    wigDict = {} #dict = {pos:{+:[+,+_e], '-':[ -,-_e]}}
    pairDict = {}
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    #-------------open GTF-----------
    if nt == 'RNA' and extend and readsType == 'PE':
        gtf_fh = open(gtf)
    #---------------------------------
    chr = ''
    for line in fh:
        lineL = line.strip().split("\t")
        name = lineL[0]
        flag = int(lineL[1])
        if chr and chr != lineL[2]:
            if readsType == 'PE' and extend:
                extendWigDict(wigDict, gtf_fh)
            outputWigDict(wigDict)
            wigDict = ''
        chr = lineL[2]
        start = int(lineL[3]) ##sam and wig are 1-based
        cigar = lineL[5] 
        regionL = computeRegion(start,cigarP.findall(cigara),name) 
        if strand:
            xs = [i[-1] for i in lineL[11:] if i.startswith('XS:A:')]
        else:
            xs = '+'
        if readsType == 'SE' and strand:
            pass
        elif readsType == 'SE' and not strand:
            pass
        elif readsType == 'PE':
            if flag & 0x2 == 2: #properly paired
                if name not in pairDict:
                    pairDict[name] = [[chr,regionL,xs]]
                else:
                    pairDict[name].append([chr,regionL,xs])
                    computeWigDict(wigDict, pairDict[name])
                    pairDict.pop(name)
                #------------------------------
            elif flag & 0x2 == 0: #unproperly paired
                for posL in regionL:
                    for pos in range(posL[0], posL[1]):
                        if pos not in wigDict:
                            wigDict[pos][xs] = [1,0]
                        else:
                            if xs not in wigDict[pos]:
                                wigDict[pos][xs] = [1,0]
                            else:
                                wigDict[pos][xs][0] += 1
                    #----------finish on region-----------
                #-----------finish all regions-------
            #--------END unproperly paired---------
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    if nt == 'RNA' and extend and readsType == 'PE':
        gtf_fh.close()
    #----last chromosome-------------------------------
    if wigDict:
        if readsType == 'PE' and extend:
            extendWigDict(wigDict, gtf_fh)
        outputWigDict(wigDict)
    #--------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


