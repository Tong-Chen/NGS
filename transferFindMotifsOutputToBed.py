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

This transfers the output of findMotifsGenome.pl or findMotifs.pl to normal bed file.
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
    usages = "%prog -i homer_out -b bed_used_for_findMotif -t RNA"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="This is the output of \
findMotifsGenome.pl or findMotifs.pl when using -find. ")
    parser.add_option("-b", "--bed-file", dest="bedin",
        metavar="inputBed", help="The original bed file which has been given to \
findMotifsGenome.pl and used to get the file given to -i. If -t is \
RNA, an at least six-column file is needed. If -t is DNA,  an at-leat \
four-column file is needed. Other columns will be ignored. \
If you run findMotifsGenome.pl using -size 200 or other parameters \
rather than <-size given>,  you may need the real bed which \
findMotifsGenome.pl used. ")
    parser.add_option("-t", "--type", dest="type",
        default='RNA', help="The attribute of given bed file \
represented sequences. For <DNA>, the strand information of regions in \
output bed file is accordant with output file of findMotifsGenome.pl. \
For <RNA>, the strand information of one region is accordant with its \
strand information in file given to -b.")
    parser.add_option("-H", "--header", dest="header",
        default=1, help="The number of header lines to skip.")
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
    bedin = options.bedin
    head = int(options.header)
    type = options.type
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    #--------get the coordinate of original bed file------
    bedinD = {}
    '''
    bedinD = {name:(chr, start, end, strand)}
    '''
    for line in open(bedin):
        lineL = line.strip().split('\t')
        key = lineL[3]
        if type == 'RNA':
            assert len(lineL) >=6, line
            if key not in bedinD:
                bedinD[key] = (lineL[0], int(lineL[1]), int(lineL[2]), lineL[5])
            else:
                print >>sys.stderr, "Diplicate region name %s" % key
        elif type == 'DNA':
            if key not in bedinD:
                bedinD[key] = (lineL[0], int(lineL[1]), int(lineL[2]))
            else:
                print >>sys.stderr, "Diplicate region name %s" % key
    #--------read findMotifGenoms.pl output---------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    motifD ={}
    '''
    motifD = {name:[(start, end, value, strand),]}
    '''
    for line in fh:
        if head:
            head -= 1
            continue
        #-----------------------------
        lineL = line.strip().split('\t')
        name = lineL[0]
        start = int(lineL[1])
        end  = start + len(lineL[2])
        motif = lineL[2]
        motif_name = lineL[3]
        if type == 'RNA':
            tmpTuple = (start,end,lineL[5], '', motif, motif_name)
        elif type == 'DNA':
            strand = lineL[4]
            if strand == '-':
                end   = int(lineL[1]) + 1
                start = end - len(lineL[2])
            tmpTuple = (start,end,lineL[5],strand, motif, motif_name)
        if name not in motifD:
            motifD[name] = []
        motifD[name].append(tmpTuple)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    #----------output-----------------
    '''
    bedinD = {name:(chr, start, end, strand)}
    motifD = {name:[(start, end, value, strand, motif, motif_name),]}
    '''
    for key, valueL in motifD.items():
        chr         = bedinD[key][0]
        genomeStart = bedinD[key][1]
        genomeEnd   = bedinD[key][2]
        if type == 'RNA':
            strand = bedinD[key][3]
        #---------------------------------
        label = 1
        for item in valueL:
            if type == 'RNA' and strand == '+':
                newstart = genomeStart + item[0]
                newend   = genomeStart + item[1]
            #pay attention to negative strand----------
            elif type == 'RNA' and strand == '-':
                newstart = genomeEnd - item[1]
                newend   = genomeEnd - item[0]
            score   = item[2]
            if type == 'DNA':
                strand = item[3]
                newstart = genomeStart + item[0]
                newend   = genomeStart + item[1]
            #----------------------------------
            if type == 'DNA':
                outputL = [chr, str(newstart), str(newend),
                        key+'@'+str(label)+item[5], item[2], item[3], item[4]]
            else:
                outputL = [chr, str(newstart), str(newend),
                        key+'@'+str(label)+item[5], item[2], strand, item[4]]
            print '\t'.join(outputL)
            label += 1
    #----------finish output------------
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



