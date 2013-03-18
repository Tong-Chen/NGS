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
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"

def flankTSS_TES(output, cs_dict, regionL, strand):
    end = cs_dict[output[0][0]]
    for seg in regionL:
        #--------------TSS---------------------------
        mrna_start = int(output[0][1])
        newstart = mrna_start - seg/2
        if newstart < 0: newstart = 0
        newend   = mrna_start + seg/2
        if newend > end: newend = end
        if strand == '+':
            name = output[0][3] + '.TSS' + str(seg)
        elif strand == '-':
            name = output[0][3] + '.TTS' + str(seg)
        output.append([output[0][0],str(newstart),str(newend),name, \
            output[0][4],output[0][5]])
        #------------TES---------------------------
        mrna_end = int(output[0][2])
        newstart = mrna_end - seg/2
        if newstart < 0: newstart = 0
        newend   = mrna_end + seg/2
        if newend > end: newend = end
        if strand == '+':
            name = output[0][3] + '.TTS' + str(seg)
        elif strand == '-':
            name = output[0][3] + '.TSS' + str(seg)
        output.append([output[0][0],str(newstart),str(newend),name, \
            output[0][4],output[0][5]])
    #--------------------------------------
#------------------------------------------------------


def outTSS_TES(output, cs_dict, regionL, strand):
    end = cs_dict[output[0][0]]
    oldseg = 0
    for seg in regionL:
        #--------------TSS---------------------------
        mrna_start = int(output[0][1])
        newstart = mrna_start - seg
        if newstart < 0: newstart = 0
        newend = mrna_start - oldseg
        if newend > 0:
            if strand == '+':
                name = output[0][3] + '.TSS-UP' + str(oldseg) + '-' + str(seg)
            elif strand == '-':
                name = output[0][3] + '.TTS-DW' + str(oldseg) + '-' + str(seg)
            output.append([output[0][0],str(newstart),str(newend),name, \
                output[0][4],output[0][5]])
        #------------TES---------------------------
        mrna_end = int(output[0][2])
        newstart = mrna_end + oldseg
        newend   = mrna_end + seg
        if newend > end: newend = end
        if newstart < end:
            if strand == '+':
                name = output[0][3] + '.TTS-DW' + str(oldseg) + '-' + str(seg) 
            elif strand == '-':
                name = output[0][3] + '.TSS-UP' + str(oldseg) + '-' + str(seg) 
            output.append([output[0][0],str(newstart),str(newend),name, \
                output[0][4],output[0][5]])
        #--------update oldseg-------------------------
        oldseg = seg
    #--------------------------------------
#------------------------------------------------------



def main():
    lensysargv = len(sys.argv)
    if lensysargv != 3:
        print >>sys.stderr, "This transfers bed 12 to bed 6 of a gene \
model. Other than the tools from Bedtools, this add name for each \
exon, intron, intergenic regions and add some TSS and TES flanking \
regions. Print the result to screen"
        print >>sys.stderr, 'Using python %s filename chrom_size' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    regionL = [500,1000,2000,5000]
    bed12 = sys.argv[1]
    chrom_size = sys.argv[2]
    cs_dict = {}
    for line in open(chrom_size):
        lineL = line.split()
        cs_dict[lineL[0]] = int(lineL[1])
    #----------------------------------------------
    name_label = 0
    for line in open(bed12):
        name_label += 1
        lineL = line.split()
        output = []
        chr    = lineL[0]
        start  = int(lineL[1])
        end    = int(lineL[2])
        name   = '_'.join([lineL[3],str(name_label)])
        value  = lineL[4]
        strand = lineL[5]
        output.append([chr, lineL[1], lineL[2], name, value, strand])
        flankTSS_TES(output, cs_dict, regionL, strand)
        outTSS_TES(output, cs_dict, regionL, strand)
        count = int(lineL[9])
        blocksize = [int(i) for i in lineL[10].split(',')[:-1]]
        blockStart = [int(i) for i in lineL[11].split(',')[:-1]]
        for i in range(count):
            if strand == '+':
                label = str(i+1)
            elif strand == '-':
                label = str(count - i)
            exon_s = start + blockStart[i]
            exon_e = exon_s + blocksize[i]
            output.append([chr,str(exon_s),str(exon_e), \
                ".".join([name,"Exon",label]),value, strand])
            intron_s = exon_e
            if i < count - 1:
                intron_e = start + blockStart[i+1]
                output.append([chr,str(intron_s),str(intron_e), \
                    ".".join([name,"Intron",label]),value, strand])
           #--------------------------------------------------------
        #-------end block---------------------------
        #if intron_s < end:
        if strand == '+':
            output.append([chr,lineL[1],lineL[6], \
                ".".join([name,"UTR5"]),value, strand])
            output.append([chr,lineL[7], lineL[2], \
                ".".join([name,"UTR3"]),value, strand])
        elif strand == '-':
            output.append([chr,lineL[1],lineL[6], \
                ".".join([name,"UTR3"]),value, strand])
            output.append([chr,lineL[7], lineL[2], \
                ".".join([name,"UTR5"]),value, strand])
        #----------------------------------------
        for i in output:
            print '\t'.join(i)
#------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


