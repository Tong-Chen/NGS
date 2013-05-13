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

'''
AI: alternative intron
AE: alternative exon
NAE: single transcript
NAI: single transcript
#only length of element larger than 2*length of boundary
LB: Left boundary
RB: Right boundary 
'''

if False:
    print "This program does not work under python 3, \
run in python 2.x."

import sys
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"

def outputBoundary(lineL, boundary):
    start = int(lineL[1])
    end = int(lineL[2])
    old = lineL[3]
    strand = lineL[5]
    if end-start > boundary * 2:
        #-------left boundary-----------
        lineL[2] = str(start+boundary)
        if strand == '+':
            lineL[3] = old + '-LB'
        elif strand == '-':
            lineL[3] = old + '-RB'
        else:
            print >>sys.stderr, \
                "Wrong strand %s " % transcriptName
            sys.exit(1)
        #---------------------------------------
        print '\t'.join(lineL)
        #------right boundary----------
        lineL[1] = str(end - boundary)
        lineL[2] = str(end)
        if strand == '+':
            lineL[3] = old + '-RB'
        elif strand == '-':
            lineL[3] = old + '-LB'
        print '\t'.join(lineL)
    #---------------end boundary------------
#------------------END------o---------------

def parse(geneDict, transcriptDict, boundary):
    for gene, transcriptL in geneDict.items():
        lentranscriptL = len(transcriptL)
        if lentranscriptL == 1:
            transcriptName = transcriptL[0]
            for type, lineL in transcriptDict[transcriptName].items():
                bed_name = lineL[3]
                if bed_name.find('Coding_exon') != -1:
                    lineL[3] = bed_name + '-NAE'
                elif bed_name.find('Intron') != -1:
                    lineL[3] = bed_name + '-NAI'
                else:
                    print >>sys.stderr,  "Wrong type %s" % bed_name
                    sys.exit(1)
                lineL[4] = gene
                print '\t'.join(lineL)
                outputBoundary(lineL,boundary)
            #--------------End ouput noAS------------
        #----------------no AS-------------------------------------
        elif lentranscriptL > 1:
            intronDict = {}
            exonDict = {}
            for tran in transcriptL:
                for key, lineL in transcriptDict[tran].items():
                    if key.find('Intron') != -1:
                        pos_key = '.'.join(lineL[:3])
                        if pos_key in intronDict:
                            intronDict[pos_key][0] += 1
                        else:
                            intronDict[pos_key] = [1]
                            intronDict[pos_key].append(lineL[4])
                    #-------------Intron------------------
                    elif key.find('Coding_exon') != -1:
                        pos_key = '.'.join(lineL[:3])
                        if pos_key in exonDict:
                            exonDict[pos_key][0] += 1
                        else:
                            exonDict[pos_key] = [1]
                            exonDict[pos_key].append(lineL[4])
                    #----------EXON-----------------
                #-----------END iterate transcriptDict--------
            #---------------END iterate all transcript---
            #-------------Begin processing intron_pos, exonDict--
            for tran in transcriptL:
                for key, lineL in transcriptDict[tran].items():
                    lineL[4] = gene
                    if key.find('Intron') != -1:
                        pos_key = '.'.join(lineL[:3])
                        time = intronDict[pos_key][0]
                        if time == lentranscriptL:
                            lineL[3] += '-'.join(['-CI', \
                                str(lentranscriptL),str(time)]) 
                            print '\t'.join(lineL)
                            outputBoundary(lineL,boundary)
                        else:
                            lineL[3] += '-'.join(['-AI', \
                                str(lentranscriptL),str(time)]) 
                            print '\t'.join(lineL)
                            outputBoundary(lineL,boundary)
                    elif key.find('Coding_exon') != -1:
                        pos_key = '.'.join(lineL[:3])
                        time = exonDict[pos_key][0]
                        if time == lentranscriptL:
                            lineL[3] += '-'.join(['-CE', \
                                str(lentranscriptL),str(time)]) 
                            print '\t'.join(lineL)
                            outputBoundary(lineL, boundary)
                        else:
                            lineL[3] += '-'.join(['-AE', \
                                str(lentranscriptL),str(time)]) 
                            print '\t'.join(lineL)
                            outputBoundary(lineL, boundary)
                    #----------------------END type-------------
                #--------END transcript iteration------------------
            #----END transcriptL-----------------------------------
        #----END all inner loop---------------------------------
    #--------------------END iteration----------------------------
#--------------------end parse--------------------


def main():
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename[- means \
sys.stdin. Usually the output of parseGTF.py.] boundary_size[default 100]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file = sys.argv[1]
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    if lensysargv > 2:
        size = int(sys.argv[2])
    else:
        size = 100
    #------------------------------------
    geneDict = {}
    transcriptDict = {}
    for line in fh:
        lineL = line.split()
        if lineL[4] == '0':
            bed_name = lineL[3]
            intron_pos = bed_name.find('.Intron')
            exon_pos = bed_name.find('.Coding_exon')
            if intron_pos != -1:
                tr_name = bed_name[:intron_pos]
                type = bed_name[intron_pos:]
            elif exon_pos != -1:
                tr_name = bed_name[:exon_pos]
                type = bed_name[exon_pos:]
            else:
                continue
            #------------------------------------
            if tr_name not in transcriptDict:
                transcriptDict[tr_name] = {}
            #-----------------------------------
            assert type not in transcriptDict[tr_name], tr_name
            transcriptDict[tr_name][type] = lineL    
        elif lineL[4] != '0':
            gene_name = '@'.join([lineL[4], lineL[0]])
            if gene_name not in geneDict:
                geneDict[gene_name] = []
            #----------------------------------
            geneDict[gene_name].append(lineL[3])
        #-----------------------------------------
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    parse(geneDict, transcriptDict, size)
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


