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
This is written to transfer BAM file to wig file. 
Currently this script only focus on Pair-end RNA-Seq reads.
One feature of this program is that it can recover the coverage for
<insert regions> between two paired reads.

For example,  we have a properly mapped reads showed below


==============================================  Genome
^^^^^^                                  $$$$$$  Two reads 
      ++++++++++++++++++++++++++++++++++        Insert regions

Traditionally, the output wig file for regions labeld '+' would be 0.
Actually,  these regions are also covered by reads. Thus,  using this
program,  you can get <true> coverage for <insert regions>.

It was also planned to transfer BAM with single-end reads to wig.
However, this is not finished yet. But other substitutaion tools or 
combined tools are available in this directory can deal with this type
of transferation.  
'''


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
expressed transcripts is preferred). )
    parser.add_option("-n", "--nucleotide-type", dest="nt",
        metavar="RNA/DNA", default="RNA", 
        help="DNA means ChIP-Seq, RNA means RNA-Seq")
    parser.add_option("-t", "--seq-type", dest="seq_Type",
        metavar="PAIREND/SINGLEEND", default="PE", 
        help="PE for pair-end reads and SE for single-end reads.")
    parser.add_option("-l", "--library-type", dest="lt",
        help="fr-unstranded,fr-firststrand,fr-secondstrand")
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

def computeWigDict(wigDict, pairL, lt):
    #wigDict = {pos:{+:[+,+_e], '-':[ -,-_e]}}
    #pairL   = [[chr,flag,[[start,end],...],xs], [chr,flag,[[start,end],...],xs]]
    #---Get relative position for two reads--------------
    # -1 means reads with smaller coordinate (more left  in Genome)
    #  0 means reads with bigger  coordinate (more right in Genome)
    # These two numbers are also used for sorting two reads and
    # extract neighbor mapped regions of two reads.
    for i in pairL:
        if i[1] & 0x10 == 0: #mate reverse
            i[1] = -1 if lt == 'fr-secondstrand' else 0
        elif i[1] & 0x20 == 0: #self reverse
            i[1] = -1 if lt == 'fr-secondstrand' else 0
    #----------------------------------------------------------- 
    assert pairL[0][0] == pairL[1][0], pairL
    assert pairL[0][3] == pairL[1][3], pairL
    xs = pairL[0][3]
    pairL.sort(key=lambda x:x[1])
    leftReadsMaxCor = pairL[0][2][-1]
    rightReadsMinCol = pairL[1][2][0]
    #--------get covergae for intervals which actually have ---
    #--coverage but no reads covering---------------------------
    #only executing when real mate inner dist larger than 0
    for pos in range(leftReadsMaxCor[1],rightReadsMinCor[0]):
        if pos not in wigDict:
            wigDict[pos][xs] = [0,-1]
        else:
            if xs not in wigDict[pos]:
                wigDict[pos][xs] = [0,-1]
            else:
                wigDict[pos][xs][1] -= 1
        #----------------END adding wigDict-----------------
    #-------END coverage for interval regions---------------
    #if leftReadsMaxCor[1] < rightReadsMinCor[0]:
    #-first get coverage for really mapped regions----
    uniqMappedPosL = [] #urgly unefficient soluable methods
    for reads in pairL:
        for region in reads[2]:
            for pos in range(region[0], region[1]):
                if pos not in uniqMappedPosL:
                    uniqMappedPosL.append(pos)
                else:
                    continue #ignore enlarge coverage for position
                             #sequences twice by one frag
                if pos not in wigDict:
                    wigDict[pos][xs] = [1,0]
                else:
                    if xs not in wigDict[pos]:
                        wigDict[pos][xs] = [1,0]
                    else:
                        wigDict[pos][xs][0] += 1
                #------------end adding to wigDict---------
            #----------------end one mapped fragments---
        #----------------end getting each mapped fragments---
    #-----------------end mapped coverage of tow reads--------------------------
#----------END computeWigDict---------------
def extendWigDict(wigDict,exonDict):
    '''Give coverage to interval regions based on the following two
    conditions.
    1.The interval located at an expressed exons.
    2.The interval has covered by other reads. If this region can not
    be sequenced randomly, it is no need to add coverage for it. If
    added,  this is no much change. 
    GTF: 1-based numbering, both closed'''
    #---get exon regions of one chromosome----------
    
#--------NED extendWigDict---------------------

def readExonRegFromGTF(gtf, lt)
    '''
    Get exon regions from GTF for each chromosome. It is such a greedy
    process that if one NT is exon in one transcript, it will be
    considered as exon [This does not always mean it will have
    coverage if interval value less than 0].

    GTF: 1-based numbering, both closed
    '''
    exonDict = {} #{chr:[[exon_s,exon_e], [], ...]} 1-based both
                    #closed
    for line in open(gtf):
        lineL = line.split('\t', 7)
        if lineL[2] == 'exon':
            xs = lineL[6] if lt != 'fr-unstranded' else '+'
            chr = lineL[0]
            if chr not in exonDict:
                exonDict[chr] = {}
            if xs not in exonDict[chr]:
                exonDict[chr][xs] = []
            exonDict[chr][xs].append([int(lineL)])
        #-----------end exon-----------
    #--------End reading file------------
    return exonDict
#-----------END read exons from GTF----------
def main():
    options = {}
    args = []
    cmdparameter(sys.argv[1:], options, args)
    #-----------------------------------
    cigarP = re.compile('([0-9]+)([A-Z])')
    file = options.filein
    output = options.fileout
    readsType = opions.seq_Type
    lt= options.lt
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
        exonDict = readExonRegFromGTF(gtf,lt)   
    #---------------------------------
    chr = ''
    for line in fh:
        lineL = line.strip().split("\t")
        name = lineL[0]
        flag = int(lineL[1])
        if chr and chr != lineL[2]:
            if readsType == 'PE' and extend:
                extendWigDict(wigDict, exonDict)
            outputWigDict(wigDict)
            wigDict = ''
        chr = lineL[2]
        start = int(lineL[3]) ##sam and wig are 1-based
        cigar = lineL[5] 
        regionL = computeRegion(start,cigarP.findall(cigara),name) 
        if lt != 'fr-unstranded':
            xs = [i[-1] for i in lineL[11:] if i.startswith('XS:A:')]
        else:
            xs = '+'
        if readsType == 'SE' and lt != 'fr-unstranded':
            pass
        elif readsType == 'SE' and lt == 'fr-unstranded': 
            pass
        elif readsType == 'PE':
            if flag & 0x2 == 2: #properly paired
                if name not in pairDict:
                    pairDict[name] = [[chr,flag,regionL,xs]]
                else:
                    pairDict[name].append([chr,flag,regionL,xs])
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
    #----last chromosome-------------------------------
    if wigDict:
        if readsType == 'PE' and extend:
            extendWigDict(wigDict, exonDict)
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


