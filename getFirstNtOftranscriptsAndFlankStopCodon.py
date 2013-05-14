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

This is designed to extract two parts of mRNA. The first part if the
first n nt regions of transcripts. [Transcriptionaal Start Site, TSS+n]
The second part is the flanking n nt regions of Stopcodon.
[Translating Stop site-n/2, Translating Stop Site+n/2].

Input file is the output of parseGTF.py.

Necessary lines:
chr12   4824595 4824625 NM_025323_15.Coding_exon.1      0       +
chr12   4827399 4827518 NM_025323_15.Coding_exon.2      0       +
chr12   4833563 4833702 NM_025323_15.Coding_exon.3      0       +
chr12   4833819 4833909 NM_025323_15.Coding_exon.4      0       +
chr12   4824413 4824595 NM_025323_15.UTR5       0       +
chr12   4833909 4834465 NM_025323_15.UTR3       0       +
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
        metavar="FILEIN", help="Ususlally input file is the output of \
parseGTF.py. Other bed files with UTR5,UTR3,Exon cooredinate in given \
format should be OK")
    parser.add_option("-t", "--tCss-length", dest="startRegion",
        metavar="200,400",
        default="200,400", help="The width you want for Transcriptional \
start site flanking regions. Single number or multiple regions separated by comma are \
accepted. Default 200,400")
    parser.add_option("-w", "--direction", dest="direction",
        default='down', metavar="down", 
        help="The ways of using <n> given to tss-length. \
[up] means extracting the upstream <n> nt of TSS. \
[down] means extracting the downstream <n> nt of TSS. \
[both] means extracting 2*<n> nt flanking TSS. Default down. \
Multiple values separated by ',' are also accepted.")
    parser.add_option("-e", "--tRts-length", dest="endRegion",
        metavar="200,400",
        default="200,400", help="The width you want for Translation \
terminating sites flanking regions. Single number or multiple regions separated by \
comma are accepted. Default 200,400.")
    parser.add_option("-o", "--orientation", dest="orientation",
        metavar='both',
        default='both', help="The ways of using <n> given to tRts-length. \
[up] means extracting the upstream <n> nt of TTS. \
[down] means extracting the downstream <n> nt of TTS. \
[both] means extracting 2*<n> nt flanking TTS. Default both. \
Multiple values separated by ',' are also accepted.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------



def processing(aDict, name, tCss_len, tCss_w, tRts_len, tRts_w):
    strand = aDict['strand']
    exonKeyL = aDict['Coding_exon'].keys()
    exonKeyL.sort()
    if strand == '+':
        if 'UTR5' not in aDict:
            for i in tCss_len:
                remainI = i
                for index in exonKeyL:
                    tmpExon  = aDict['Coding_exon'][index]
                    in_start = int(tmpExon[1])
                    in_end   = int(tmpExon[2])
                    in_width = in_end - in_start
                    if remainI > in_width:
                        tmpL = [tmpExon[0], str(in_start), str(in_end),
                                name+'.TcSS.dw'+str(i),tmpExon[4],
                                strand]
                        print '\t'.join(tmpL)
                        remainI = remainI - in_width
                    else:
                        tmpL = [tmpExon[0], str(in_start),
                                str(in_start+remainI),
                                name+'.TcSS.dw'+str(i),tmpExon[4],
                                strand]
                        print '\t'.join(tmpL)
                        remainI = 0
                        break
                #-------------Get full length-----------------
            #-----------------END one type of length---------------------
        #-------END no UTR5----------------------------------------------
        else:
            for i in tCss_len:
                remainI = i
                #print aDict['UTR5']
                for UTR5 in aDict['UTR5']:
                    #print UTR5
                    start = int(UTR5[1])
                    end   = int(UTR5[2])
                    width = end -start
                    if remainI <= width:
                        tmpL = [UTR5[0], str(start), str(start+remainI),
                                name+'.TcSS.dw'+str(i),UTR5[4],
                                strand]
                        print '\t'.join(tmpL)
                        remainI = 0
                        break
                    else:
                        tmpL = [UTR5[0], str(start), str(end),
                                name+'.TcSS.dw'+str(i),UTR5[4],
                                strand]
                        print '\t'.join(tmpL)
                        remainI = remainI - width 
                #---------Get seg length from UTR5-----------
                #---------Get seg length from Exon-----------
                if remainI == 0:
                    continue
                for index in exonKeyL:
                    tmpExon  = aDict['Coding_exon'][index]
                    in_start = int(tmpExon[1])
                    in_end   = int(tmpExon[2])
                    in_width = in_end - in_start
                    if remainI > in_width:
                        tmpL = [tmpExon[0], str(in_start), str(in_end),
                                name+'.TcSS.dw'+str(i),tmpExon[4],
                                strand]
                        print '\t'.join(tmpL)
                        remainI = remainI - in_width
                    else:
                        tmpL = [tmpExon[0], str(in_start),
                                str(in_start+remainI),
                                name+'.TcSS.dw'+str(i),tmpExon[4],
                                strand]
                        print '\t'.join(tmpL)
                        remainI = 0
                        break
                #-------------Get full length-----------------
            #-----------------END one type of length---------------------
        #-------END no UTR5----------------------------------------------
        #--------------------BEGIN UTR3---------------------------------
        exonKeyL.reverse()
        if 'UTR3' not in aDict:
            for i in tRts_len:
                remainI = i
                for type in tRts_w:
                    if type == 'down' :
                        continue
                    for index in exonKeyL:
                        tmpExon  = aDict['Coding_exon'][index]
                        in_start = int(tmpExon[1])
                        in_end   = int(tmpExon[2])
                        in_width = in_end - in_start
                        if remainI > in_width:
                            tmpL = [tmpExon[0], str(in_start), str(in_end),
                                    name+'.TsTS.'+type+str(i),tmpExon[4],
                                    strand]
                            print '\t'.join(tmpL)
                            remainI = remainI - in_width
                        else:
                            tmpL = [tmpExon[0], str(in_end-remainI),
                                    str(in_end),
                                    name+'.TsTS.'+type+str(i),tmpExon[4],
                                    strand]
                            print '\t'.join(tmpL)
                            remainI = 0
                            break
                    #-------------Get full length-----------------
                #-----------------END one type of length---------------------
            #-----------------END various length---------------------
        #-------END no UTR3----------------------------------------------
        else:
            for i in tRts_len:
                for type in tRts_w:
                    remainI = i
                    if type == 'UP':
                        for index in exonKeyL:
                            tmpExon  = aDict['Coding_exon'][index]
                            in_start = int(tmpExon[1])
                            in_end   = int(tmpExon[2])
                            in_width = in_end - in_start
                            if remainI > in_width:
                                tmpL = [tmpExon[0], str(in_start), str(in_end),
                                        name+'.TsTS.'+type+str(i),tmpExon[4],
                                        strand]
                                print '\t'.join(tmpL)
                                remainI = remainI - in_width
                            else:
                                tmpL = [tmpExon[0], str(in_end-remainI),
                                        str(in_end),
                                        name+'.TsTS.'+type+str(i),tmpExon[4],
                                        strand]
                                print '\t'.join(tmpL)
                                remainI = 0
                                break
                        #-------------Get full length-----------------
                    elif type == 'DW':
                        for UTR3 in aDict['UTR3']:
                            in_start = int(UTR3[1])
                            in_end   = int(UTR3[2])
                            in_width = in_end - in_start
                            if remainI > in_width:
                                tmpL = [UTR3[0],UTR3[1],UTR3[2], 
                                    name+'.TsTS.DW'+str(i),UTR3[4],strand]
                                remainI = remainI - in_width
                                print '\t'.join(tmpL)
                            else:
                                tmpL = [UTR3[0],UTR3[1],str(in_start+remainI), 
                                    name+'.TsTS.DW'+str(i),UTR3[4],strand]
                                print '\t'.join(tmpL)
                                remainI = 0
                                break
                            #----------------------------------------------
                        #------------------------------------
                    elif type == 'both':
                        exonRemainI = remainI
                        UTR3RemainI = remainI
                        for index in exonKeyL:
                            tmpExon  = aDict['Coding_exon'][index]
                            in_start = int(tmpExon[1])
                            in_end   = int(tmpExon[2])
                            in_width = in_end - in_start
                            if exonRemainI > in_width:
                                tmpL = [tmpExon[0], str(in_start), str(in_end),
                                        name+'.TsTS.'+type+str(i),tmpExon[4],
                                        strand]
                                print '\t'.join(tmpL)
                                exonRemainI = exonRemainI - in_width
                            else:
                                tmpL = [tmpExon[0], str(in_end-exonRemainI),
                                        str(in_end),
                                        name+'.TsTS.'+type+str(i),tmpExon[4],
                                        strand]
                                print '\t'.join(tmpL)
                                exonRemainI = 0
                                break
                        #-------------Get full length-----------------
                        for UTR3 in aDict['UTR3']:
                            in_start = int(UTR3[1])
                            in_end   = int(UTR3[2])
                            in_width = in_end - in_start
                            if UTR3RemainI > in_width:
                                tmpL = [UTR3[0],UTR3[1],UTR3[2], 
                                    name+'.TsTS.'+type+str(i),UTR3[4],strand]
                                UTR3RemainI = UTR3RemainI - in_width
                                print '\t'.join(tmpL)
                            else:
                                tmpL = \
                                    [UTR3[0],UTR3[1],str(in_start+UTR3RemainI), 
                                    name+'.TsTS.'+type+str(i),UTR3[4],strand]
                                print '\t'.join(tmpL)
                                UTR3RemainI = 0
                                break
                            #----------------------------------------------
                        #------------------------------------
                #-----------------END one type of length---------------------
            #-----------------END various length---------------------
        #--------------------END UTR3---------------------
    #----------------------END all positive strand------
    elif strand == '-':
        exonKeyL.reverse()
        if 'UTR5' not in aDict:
            for i in tCss_len:
                remainI = i
                for index in exonKeyL:
                    tmpExon  = aDict['Coding_exon'][index]
                    in_start = int(tmpExon[1])
                    in_end   = int(tmpExon[2])
                    in_width = in_end - in_start
                    if remainI > in_width:
                        tmpL = [tmpExon[0], str(in_start), str(in_end),
                                name+'.TcSS.dw'+str(i),tmpExon[4],
                                strand]
                        print '\t'.join(tmpL)
                        remainI = remainI - in_width
                    else:
                        tmpL = [tmpExon[0], str(in_end-remainI),
                                str(in_end),
                                name+'.TcSS.dw'+str(i),tmpExon[4],
                                strand]
                        print '\t'.join(tmpL)
                        remainI = 0
                        break
                #-------------Get full length-----------------
            #-----------------END one type of length---------------------
        #-------END no UTR5----------------------------------------------
        else:
            for i in tCss_len:
                remainI = i
                UTR5L = aDict['UTR5']
                UTR5L.reverse()
                for UTR5 in UTR5L:
                    start = int(UTR5[1])
                    end   = int(UTR5[2])
                    width = end -start
                    if remainI <= width:
                        tmpL = [UTR5[0], str(end-remainI), str(end),
                                name+'.TcSS.dw'+str(i),UTR5[4],
                                strand]
                        print '\t'.join(tmpL)
                        remainI = 0
                        break
                    else:
                        tmpL = [UTR5[0], str(start), str(end),
                                name+'.TcSS.dw'+str(i),UTR5[4],
                                strand]
                        print '\t'.join(tmpL)
                        remainI = remainI - width 
                #---------Get seg length from UTR5-----------
                #---------Get seg length from Exon-----------
                if remainI == 0:
                    continue
                for index in exonKeyL:
                    tmpExon  = aDict['Coding_exon'][index]
                    in_start = int(tmpExon[1])
                    in_end   = int(tmpExon[2])
                    in_width = in_end - in_start
                    if remainI > in_width:
                        tmpL = [tmpExon[0], str(in_start), str(in_end),
                                name+'.TcSS.dw'+str(i),tmpExon[4],
                                strand]
                        print '\t'.join(tmpL)
                        remainI = remainI - in_width
                    else:
                        tmpL = [tmpExon[0], str(in_end-remainI),
                                str(in_end),
                                name+'.TcSS.dw'+str(i),tmpExon[4],
                                strand]
                        print '\t'.join(tmpL)
                        remainI = 0
                        break
                    #-------------Get full length-----------------
                #-------------Get full length-----------------
            #-----------------END one type of length---------------------
        #-------END no UTR5----------------------------------------------
        #--------------------BEGIN UTR3---------------------------------
        exonKeyL.reverse()
        if "UTR3" not in aDict:
            for i in tRts_len:
                remainI = i
                for type in tRts_w:
                    if type == 'down' :
                        continue
                    for index in exonKeyL:
                        tmpExon  = aDict['Coding_exon'][index]
                        in_start = int(tmpExon[1])
                        in_end   = int(tmpExon[2])
                        in_width = in_end - in_start
                        if remainI > in_width:
                            tmpL = [tmpExon[0], str(in_start), str(in_end),
                                    name+'.TsTS.'+type+str(i),tmpExon[4],
                                    strand]
                            print '\t'.join(tmpL)
                            remainI = remainI - in_width
                        else:
                            tmpL = [tmpExon[0], str(in_start),
                                    str(in_start+remainI),
                                    name+'.TsTS.'+type+str(i),tmpExon[4],
                                    strand]
                            print '\t'.join(tmpL)
                            remainI = 0
                            break
                    #-------------Get full length-----------------
                #-----------------END one type of length---------------------
            #-----------------END various length---------------------
        #-------END no UTR3----------------------------------------------
        else:
            UTR3L = aDict['UTR3']
            UTR3L.reverse()
            for i in tRts_len:
                for type in tRts_w:
                    remainI = i
                    if type == 'UP':
                        for index in exonKeyL:
                            tmpExon  = aDict['Coding_exon'][index]
                            in_start = int(tmpExon[1])
                            in_end   = int(tmpExon[2])
                            in_width = in_end - in_start
                            if remainI > in_width:
                                tmpL = [tmpExon[0], str(in_start), str(in_end),
                                        name+'.TsTS.'+type+str(i),tmpExon[4],
                                        strand]
                                print '\t'.join(tmpL)
                                remainI = remainI - in_width
                            else:
                                tmpL = [tmpExon[0], str(in_start),
                                        str(in_start+remainI),
                                        name+'.TsTS.'+type+str(i),tmpExon[4],
                                        strand]
                                print '\t'.join(tmpL)
                                remainI = 0
                                break
                        #-------------Get full length-----------------
                    elif type == 'DW':
                        for UTR3 in UTR3L:
                            in_start = int(UTR3[1])
                            in_end   = int(UTR3[2])
                            in_width = in_end - in_start
                            if remainI > in_width:
                                tmpL = [UTR3[0],UTR3[1],UTR3[2], 
                                    name+'.TsTS.DW'+str(i),UTR3[4],strand]
                                remainI = remainI - in_width
                                print '\t'.join(tmpL)
                            else:
                                tmpL = [UTR3[0],in_end-remainI,UTR3[2], 
                                    name+'.TsTS.DW'+str(i),UTR3[4],strand]
                                print '\t'.join(tmpL)
                                remainI = 0
                                break
                            #----------------------------------------------
                        #------------------------------------
                    elif type == 'both':
                        exonRemainI = remainI
                        UTR3RemainI = remainI
                        for index in exonKeyL:
                            tmpExon  = aDict['Coding_exon'][index]
                            in_start = int(tmpExon[1])
                            in_end   = int(tmpExon[2])
                            in_width = in_end - in_start
                            if exonRemainI > in_width:
                                tmpL = [tmpExon[0], str(in_start), str(in_end),
                                        name+'.TsTS.'+type+str(i),tmpExon[4],
                                        strand]
                                print '\t'.join(tmpL)
                                exonRemainI = exonRemainI - in_width
                            else:
                                tmpL = [tmpExon[0], str(in_start),
                                        str(in_start+exonRemainI),
                                        name+'.TsTS.'+type+str(i),tmpExon[4],
                                        strand]
                                print '\t'.join(tmpL)
                                exonRemainI = 0
                                break
                        #-------------Get full length-----------------
                        for UTR3 in UTR3L:
                            in_start = int(UTR3[1])
                            in_end   = int(UTR3[2])
                            in_width = in_end - in_start
                            if UTR3RemainI > in_width:
                                tmpL = [UTR3[0],UTR3[1],UTR3[2], 
                                    name+'.TsTS.'+type+str(i),UTR3[4],strand]
                                UTR3RemainI = UTR3RemainI - in_width
                                print '\t'.join(tmpL)
                            else:
                                tmpL = \
                                    [UTR3[0],str(in_end-UTR3RemainI),UTR3[2], 
                                    name+'.TsTS.'+type+str(i),UTR3[4],strand]
                                print '\t'.join(tmpL)
                                UTR3RemainI = 0
                                break
                            #----------------------------------------------
                        #------------------------------------
                    #----------------------------------
                #-----------------END one type of length---------------------
            #-----------------END various length---------------------
        #----------------------END all positive strand------
    #---------------------------------------------
#-------------------END processing------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    tCss_len = [int(i) for i in options.startRegion.split(',')]
    tCss_w = options.direction.split(',') 
    tRts_len = [int(i) for i in options.endRegion.split(',')]
    tRts_w = options.orientation.split(',')
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    key = ''
    aDict = {}
    '''
    aDict - {'Coding_exon':{1:[chr1, start, end], 2:[chr2, start,
            end]}, 'UTR5':[[], []], 'UTR3':[[], []]}
    '''
    for line in fh:
        lineL = line.split()
        forthCol = lineL[3]
        newkey = forthCol.split('.')[0] 
        if key and newkey != key:
            processing(aDict, key, tCss_len, tCss_w, tRts_len, tRts_w)
            aDict = {}
        key = newkey
        if 'strand' not in aDict:
            aDict['strand'] = lineL[5]
        else:
            assert lineL[5] == aDict['strand']
        if forthCol.find('exon') != -1:
            null,type,num = forthCol.split('.')
            num = int(num)
            if type not in aDict:
                aDict[type] = {}
            aDict[type][num] = lineL
        elif forthCol.find('UTR5') != -1 or forthCol.find('UTR3') != -1:
            type = forthCol.split('.',1)[1] 
            if type not in aDict:
                aDict[type] = [lineL]
            else:
                aDict[type].append(lineL)
            #----------------------------------------------------------
        #------------------------------------------------------
    #-------------END reading file----------
    if key:
        processing(aDict, key, tCss_len, tCss_w, tRts_len, tRts_w)
        aDict = {}
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


#'''test.bed-----------------------------
#chr18	38410053	38410739	NM_024179_17.Intron.1	0	+
#chr18	38410022	38410053	NM_024179_17.Coding_exon.1	0	+
#chr18	38410854	38410936	NM_024179_17.Intron.2	0	+
#chr18	38411054	38412474	NM_024179_17.Intron.3	0	+
#chr18	38412619	38413561	NM_024179_17.Intron.4	0	+
#chr18	38413723	38413933	NM_024179_17.Intron.5	0	+
#chr18	38414022	38414123	NM_024179_17.Intron.6	0	+
#chr18	38414220	38416911	NM_024179_17.Intron.7	0	+
#chr18	38417054	38417681	NM_024179_17.Intron.8	0	+
#chr18	38417840	38417923	NM_024179_17.Intron.9	0	+
#chr18	38418016	38419550	NM_024179_17.Intron.10	0	+
#chr18	38419710	38420717	NM_024179_17.Intron.11	0	+
#chr18	38410739	38410854	NM_024179_17.Coding_exon.2	0	+
#chr18	38410936	38411054	NM_024179_17.Coding_exon.3	0	+
#chr18	38412474	38412619	NM_024179_17.Coding_exon.4	0	+
#chr18	38413561	38413723	NM_024179_17.Coding_exon.5	0	+
#chr18	38413933	38414022	NM_024179_17.Coding_exon.6	0	+
#chr18	38414123	38414220	NM_024179_17.Coding_exon.7	0	+
#chr18	38416911	38417054	NM_024179_17.Coding_exon.8	0	+
#chr18	38417681	38417840	NM_024179_17.Coding_exon.9	0	+
#chr18	38417923	38418016	NM_024179_17.Coding_exon.10	0	+
#chr18	38419550	38419710	NM_024179_17.Coding_exon.11	0	+
#chr18	38420717	38420938	NM_024179_17.Coding_exon.12	0	+
#chr18	38420939	38421533	NM_024179_17.Intron.12	0	+
#chr18	38409902	38422283	NM_024179_17	0610009O20Rik	+
#chr18	38409652	38410152	NM_024179_17.TSS500	0	+
#chr18	38422032	38422532	NM_024179_17.TTS500	0	+
#chr18	38409402	38410402	NM_024179_17.TSS1000	0	+
#chr18	38421782	38422782	NM_024179_17.TTS1000	0	+
#chr18	38408902	38410902	NM_024179_17.TSS2000	0	+
#chr18	38421282	38423282	NM_024179_17.TTS2000	0	+
#chr18	38407402	38412402	NM_024179_17.TSS5000	0	+
#chr18	38419782	38424782	NM_024179_17.TTS5000	0	+
#chr18	38409402	38409902	NM_024179_17.TSS-UP0-500	0	+
#chr18	38422282	38422782	NM_024179_17.TTS-DW0-500	0	+
#chr18	38408902	38409402	NM_024179_17.TSS-UP500-1000	0	+
#chr18	38422782	38423282	NM_024179_17.TTS-DW500-1000	0	+
#chr18	38407902	38408902	NM_024179_17.TSS-UP1000-2000	0	+
#chr18	38423282	38424282	NM_024179_17.TTS-DW1000-2000	0	+
#chr18	38404902	38407902	NM_024179_17.TSS-UP2000-5000	0	+
#chr18	38424282	38427282	NM_024179_17.TTS-DW2000-5000	0	+
#chr18	38409902	38410022	NM_024179_17.UTR5	0	+
#chr18	38420938	38420939	NM_024179_17.UTR3	0	+
#chr18	38421533	38422283	NM_024179_17.UTR3	0	+
#chr2	175017774	175147602	NM_001177543_18.Intron.1	0	-
#chr2	175151549	175152706	NM_001177543_18.Intron.2	0	-
#chr2	175150243	175151549	NM_001177543_18.Coding_exon.1	0	-
#chr2	175152767	175152972	NM_001177543_18.Intron.3	0	-
#chr2	175152706	175152767	NM_001177543_18.Coding_exon.2	0	-
#chr2	175152972	175153099	NM_001177543_18.Coding_exon.3	0	-
#chr2	175153099	175160123	NM_001177543_18.Intron.4	0	-
#chr2	175160123	175160126	NM_001177543_18.Coding_exon.4	0	-
#chr2	175160178	175163633	NM_001177543_18.Intron.5	0	-
#chr2	175017505	175163713	NM_001177543_18	0610010B08Rik	-
#chr2	175017255	175017755	NM_001177543_18.TTS500	0	-
#chr2	175163462	175163962	NM_001177543_18.TSS500	0	-
#chr2	175017005	175018005	NM_001177543_18.TTS1000	0	-
#chr2	175163212	175164212	NM_001177543_18.TSS1000	0	-
#chr2	175016505	175018505	NM_001177543_18.TTS2000	0	-
#chr2	175162712	175164712	NM_001177543_18.TSS2000	0	-
#chr2	175015005	175020005	NM_001177543_18.TTS5000	0	-
#chr2	175161212	175166212	NM_001177543_18.TSS5000	0	-
#chr2	175017005	175017505	NM_001177543_18.TTS-DW0-500	0	-
#chr2	175163712	175164212	NM_001177543_18.TSS-UP0-500	0	-
#chr2	175016505	175017005	NM_001177543_18.TTS-DW500-1000	0	-
#chr2	175164212	175164712	NM_001177543_18.TSS-UP500-1000	0	-
#chr2	175015505	175016505	NM_001177543_18.TTS-DW1000-2000	0	-
#chr2	175164712	175165712	NM_001177543_18.TSS-UP1000-2000	0	-
#chr2	175012505	175015505	NM_001177543_18.TTS-DW2000-5000	0	-
#chr2	175165712	175168712	NM_001177543_18.TSS-UP2000-5000	0	-
#chr2	175017505	175017774	NM_001177543_18.UTR3	0	-
#chr2	175147602	175150243	NM_001177543_18.UTR3	0	-
#chr2	175160126	175160178	NM_001177543_18.UTR5	0	-
#chr2	175163633	175163713	NM_001177543_18.UTR5	0	-
#chr11	120210176	120212070	NR_038126_16.Intron.1	0	+
#chr11	120209991	120210176	NR_038126_16.Coding_exon.1	0	+
#chr11	120212070	120212504	NR_038126_16.Coding_exon.2	0	+
#chr11	120209991	120212504	NR_038126_16	0610009L18Rik	+
#chr11	120209741	120210241	NR_038126_16.TSS500	0	+
#chr11	120212253	120212753	NR_038126_16.TTS500	0	+
#chr11	120209491	120210491	NR_038126_16.TSS1000	0	+
#chr11	120212003	120213003	NR_038126_16.TTS1000	0	+
#chr11	120208991	120210991	NR_038126_16.TSS2000	0	+
#chr11	120211503	120213503	NR_038126_16.TTS2000	0	+
#chr11	120207491	120212491	NR_038126_16.TSS5000	0	+
#chr11	120210003	120215003	NR_038126_16.TTS5000	0	+
#chr11	120209491	120209991	NR_038126_16.TSS-UP0-500	0	+
#chr11	120212503	120213003	NR_038126_16.TTS-DW0-500	0	+
#chr11	120208991	120209491	NR_038126_16.TSS-UP500-1000	0	+
#chr11	120213003	120213503	NR_038126_16.TTS-DW500-1000	0	+
#chr11	120207991	120208991	NR_038126_16.TSS-UP1000-2000	0	+
#chr11	120213503	120214503	NR_038126_16.TTS-DW1000-2000	0	+
#chr11	120204991	120207991	NR_038126_16.TSS-UP2000-5000	0	+
#chr11	120214503	120217503	NR_038126_16.TTS-DW2000-5000	0	+
#chr5	130695789	130705318	NM_001081394_4.Intron.1	0	+
#chr5	130705496	130707624	NM_001081394_4.Intron.2	0	+
#chr5	130705337	130705496	NM_001081394_4.Coding_exon.1	0	+
#chr5	130707910	130712672	NM_001081394_4.Intron.3	0	+
#chr5	130712823	130714116	NM_001081394_4.Intron.4	0	+
#chr5	130714300	130716217	NM_001081394_4.Intron.5	0	+
#chr5	130716361	130717147	NM_001081394_4.Intron.6	0	+
#chr5	130707624	130707910	NM_001081394_4.Coding_exon.2	0	+
#chr5	130712672	130712823	NM_001081394_4.Coding_exon.3	0	+
#chr5	130714116	130714300	NM_001081394_4.Coding_exon.4	0	+
#chr5	130716217	130716361	NM_001081394_4.Coding_exon.5	0	+
#chr5	130717147	130717168	NM_001081394_4.Coding_exon.6	0	+
#chr5	130695613	130719635	NM_001081394_4	0610007L01Rik	+
#chr5	130695363	130695863	NM_001081394_4.TSS500	0	+
#chr5	130719384	130719884	NM_001081394_4.TTS500	0	+
#chr5	130695113	130696113	NM_001081394_4.TSS1000	0	+
#chr5	130719134	130720134	NM_001081394_4.TTS1000	0	+
#chr5	130694613	130696613	NM_001081394_4.TSS2000	0	+
#chr5	130718634	130720634	NM_001081394_4.TTS2000	0	+
#chr5	130693113	130698113	NM_001081394_4.TSS5000	0	+
#chr5	130717134	130722134	NM_001081394_4.TTS5000	0	+
#chr5	130695113	130695613	NM_001081394_4.TSS-UP0-500	0	+
#chr5	130719634	130720134	NM_001081394_4.TTS-DW0-500	0	+
#chr5	130694613	130695113	NM_001081394_4.TSS-UP500-1000	0	+
#chr5	130720134	130720634	NM_001081394_4.TTS-DW500-1000	0	+
#chr5	130693613	130694613	NM_001081394_4.TSS-UP1000-2000	0	+
#chr5	130720634	130721634	NM_001081394_4.TTS-DW1000-2000	0	+
#chr5	130690613	130693613	NM_001081394_4.TSS-UP2000-5000	0	+
#chr5	130721634	130724634	NM_001081394_4.TTS-DW2000-5000	0	+
#chr5	130695613	130695789	NM_001081394_4.UTR5	0	+
#chr5	130705318	130705337	NM_001081394_4.UTR5	0	+
#chr5	130717168	130719635	NM_001081394_4.UTR3	0	+
#chr7	52823749	52829782	NR_038165_1.Intron.1	0	-
#chr7	52829892	52829977	NR_038165_1.Intron.2	0	-
#chr7	52830147	52830496	NR_038165_1.Intron.3	0	-
#chr7	52823164	52823749	NR_038165_1.Coding_exon.1	0	-
#chr7	52829782	52829892	NR_038165_1.Coding_exon.2	0	-
#chr7	52829977	52830147	NR_038165_1.Coding_exon.3	0	-
#chr7	52830496	52830546	NR_038165_1.Coding_exon.4	0	-
#chr7	52823164	52830546	NR_038165_1	0610005C13Rik	-
#chr7	52822914	52823414	NR_038165_1.TTS500	0	-
#chr7	52830295	52830795	NR_038165_1.TSS500	0	-
#chr7	52822664	52823664	NR_038165_1.TTS1000	0	-
#chr7	52830045	52831045	NR_038165_1.TSS1000	0	-
#chr7	52822164	52824164	NR_038165_1.TTS2000	0	-
#chr7	52829545	52831545	NR_038165_1.TSS2000	0	-
#chr7	52820664	52825664	NR_038165_1.TTS5000	0	-
#chr7	52828045	52833045	NR_038165_1.TSS5000	0	-
#chr7	52822664	52823164	NR_038165_1.TTS-DW0-500	0	-
#chr7	52830545	52831045	NR_038165_1.TSS-UP0-500	0	-
#chr7	52822164	52822664	NR_038165_1.TTS-DW500-1000	0	-
#chr7	52831045	52831545	NR_038165_1.TSS-UP500-1000	0	-
#chr7	52821164	52822164	NR_038165_1.TTS-DW1000-2000	0	-
#chr7	52831545	52832545	NR_038165_1.TSS-UP1000-2000	0	-
#chr7	52818164	52821164	NR_038165_1.TTS-DW2000-5000	0	-
#chr7	52832545	52835545	NR_038165_1.TSS-UP2000-5000	0	-
#chr7	52823749	52826355	NR_038166_2.Intron.1	0	-
#chr7	52826562	52829782	NR_038166_2.Intron.2	0	-
#chr7	52829892	52829977	NR_038166_2.Intron.3	0	-
#chr7	52830147	52830496	NR_038166_2.Intron.4	0	-
#chr7	52823164	52823749	NR_038166_2.Coding_exon.1	0	-
#chr7	52826355	52826562	NR_038166_2.Coding_exon.2	0	-
#chr7	52829782	52829892	NR_038166_2.Coding_exon.3	0	-
#chr7	52829977	52830147	NR_038166_2.Coding_exon.4	0	-
#chr7	52830496	52830546	NR_038166_2.Coding_exon.5	0	-
#chr7	52823164	52830546	NR_038166_2	0610005C13Rik	-
#chr7	52822914	52823414	NR_038166_2.TTS500	0	-
#chr7	52830295	52830795	NR_038166_2.TSS500	0	-
#chr7	52822664	52823664	NR_038166_2.TTS1000	0	-
#chr7	52830045	52831045	NR_038166_2.TSS1000	0	-
#chr7	52822164	52824164	NR_038166_2.TTS2000	0	-
#chr7	52829545	52831545	NR_038166_2.TSS2000	0	-
#chr7	52820664	52825664	NR_038166_2.TTS5000	0	-
#chr7	52828045	52833045	NR_038166_2.TSS5000	0	-
#chr7	52822664	52823164	NR_038166_2.TTS-DW0-500	0	-
#chr7	52830545	52831045	NR_038166_2.TSS-UP0-500	0	-
#chr7	52822164	52822664	NR_038166_2.TTS-DW500-1000	0	-
#chr7	52831045	52831545	NR_038166_2.TSS-UP500-1000	0	-
#chr7	52821164	52822164	NR_038166_2.TTS-DW1000-2000	0	-
#chr7	52831545	52832545	NR_038166_2.TSS-UP1000-2000	0	-
#chr7	52818164	52821164	NR_038166_2.TTS-DW2000-5000	0	-
#chr7	52832545	52835545	NR_038166_2.TSS-UP2000-5000	0	-
#chr5	31351129	31351834	NM_027855_3.Intron.1	0	+
#chr5	31351045	31351129	NM_027855_3.Coding_exon.1	0	+
#chr5	31351953	31354569	NM_027855_3.Intron.2	0	+
#chr5	31354641	31354834	NM_027855_3.Intron.3	0	+
#chr5	31354906	31355135	NM_027855_3.Intron.4	0	+
#chr5	31355257	31356333	NM_027855_3.Intron.5	0	+
#chr5	31356431	31356635	NM_027855_3.Intron.6	0	+
#chr5	31351834	31351953	NM_027855_3.Coding_exon.2	0	+
#chr5	31354569	31354641	NM_027855_3.Coding_exon.3	0	+
#chr5	31354834	31354906	NM_027855_3.Coding_exon.4	0	+
#chr5	31355135	31355257	NM_027855_3.Coding_exon.5	0	+
#chr5	31356333	31356431	NM_027855_3.Coding_exon.6	0	+
#chr5	31356635	31356740	NM_027855_3.Coding_exon.7	0	+
#chr5	31351012	31356996	NM_027855_3	0610007C21Rik	+
#chr5	31350762	31351262	NM_027855_3.TSS500	0	+
#chr5	31356745	31357245	NM_027855_3.TTS500	0	+
#chr5	31350512	31351512	NM_027855_3.TSS1000	0	+
#chr5	31356495	31357495	NM_027855_3.TTS1000	0	+
#chr5	31350012	31352012	NM_027855_3.TSS2000	0	+
#chr5	31355995	31357995	NM_027855_3.TTS2000	0	+
#chr5	31348512	31353512	NM_027855_3.TSS5000	0	+
#chr5	31354495	31359495	NM_027855_3.TTS5000	0	+
#chr5	31350512	31351012	NM_027855_3.TSS-UP0-500	0	+
#chr5	31356995	31357495	NM_027855_3.TTS-DW0-500	0	+
#chr5	31350012	31350512	NM_027855_3.TSS-UP500-1000	0	+
#chr5	31357495	31357995	NM_027855_3.TTS-DW500-1000	0	+
#chr5	31349012	31350012	NM_027855_3.TSS-UP1000-2000	0	+
#chr5	31357995	31358995	NM_027855_3.TTS-DW1000-2000	0	+
#chr5	31346012	31349012	NM_027855_3.TSS-UP2000-5000	0	+
#chr5	31358995	31361995	NM_027855_3.TTS-DW2000-5000	0	+
#chr5	31351012	31351045	NM_027855_3.UTR5	0	+
#chr5	31356740	31356996	NM_027855_3.UTR3	0	+
#chr11	51499148	51499571	NM_025319_14.Coding_exon.1	0	-
#chr11	51499594	51502049	NM_025319_14.Intron.1	0	-
#chr11	51498886	51502136	NM_025319_14	0610009B22Rik	-
#chr11	51498636	51499136	NM_025319_14.TTS500	0	-
#chr11	51501885	51502385	NM_025319_14.TSS500	0	-
#chr11	51498386	51499386	NM_025319_14.TTS1000	0	-
#chr11	51501635	51502635	NM_025319_14.TSS1000	0	-
#chr11	51497886	51499886	NM_025319_14.TTS2000	0	-
#chr11	51501135	51503135	NM_025319_14.TSS2000	0	-
#chr11	51496386	51501386	NM_025319_14.TTS5000	0	-
#chr11	51499635	51504635	NM_025319_14.TSS5000	0	-
#chr11	51498386	51498886	NM_025319_14.TTS-DW0-500	0	-
#chr11	51502135	51502635	NM_025319_14.TSS-UP0-500	0	-
#chr11	51497886	51498386	NM_025319_14.TTS-DW500-1000	0	-
#chr11	51502635	51503135	NM_025319_14.TSS-UP500-1000	0	-
#chr11	51496886	51497886	NM_025319_14.TTS-DW1000-2000	0	-
#chr11	51503135	51504135	NM_025319_14.TSS-UP1000-2000	0	-
#chr11	51493886	51496886	NM_025319_14.TTS-DW2000-5000	0	-
#chr11	51504135	51507135	NM_025319_14.TSS-UP2000-5000	0	-
#chr11	51498886	51499148	NM_025319_14.UTR3	0	-
#chr11	51499571	51499594	NM_025319_14.UTR5	0	-
#chr11	51502049	51502136	NM_025319_14.UTR5	0	-
#chr2	176266848	176270305	NM_001177543_5_22.Intron.1	0	+
#chr2	176270360	176277384	NM_001177543_5_22.Intron.2	0	+
#chr2	176270357	176270360	NM_001177543_5_22.Coding_exon.1	0	+
#chr2	176277511	176277716	NM_001177543_5_22.Intron.3	0	+
#chr2	176277777	176278934	NM_001177543_5_22.Intron.4	0	+
#chr2	176277384	176277511	NM_001177543_5_22.Coding_exon.2	0	+
#chr2	176277716	176277777	NM_001177543_5_22.Coding_exon.3	0	+
#chr2	176278934	176280240	NM_001177543_5_22.Coding_exon.4	0	+
#chr2	176282881	176598568	NM_001177543_5_22.Intron.5	0	+
#chr2	176266768	176598837	NM_001177543_5_22	0610010B08Rik	+
#chr2	176266518	176267018	NM_001177543_5_22.TSS500	0	+
#chr2	176598586	176599086	NM_001177543_5_22.TTS500	0	+
#chr2	176266268	176267268	NM_001177543_5_22.TSS1000	0	+
#chr2	176598336	176599336	NM_001177543_5_22.TTS1000	0	+
#chr2	176265768	176267768	NM_001177543_5_22.TSS2000	0	+
#chr2	176597836	176599836	NM_001177543_5_22.TTS2000	0	+
#chr2	176264268	176269268	NM_001177543_5_22.TSS5000	0	+
#chr2	176596336	176601336	NM_001177543_5_22.TTS5000	0	+
#chr2	176266268	176266768	NM_001177543_5_22.TSS-UP0-500	0	+
#chr2	176598836	176599336	NM_001177543_5_22.TTS-DW0-500	0	+
#chr2	176265768	176266268	NM_001177543_5_22.TSS-UP500-1000	0	+
#chr2	176599336	176599836	NM_001177543_5_22.TTS-DW500-1000	0	+
#chr2	176264768	176265768	NM_001177543_5_22.TSS-UP1000-2000	0	+
#chr2	176599836	176600836	NM_001177543_5_22.TTS-DW1000-2000	0	+
#chr2	176261768	176264768	NM_001177543_5_22.TSS-UP2000-5000	0	+
#chr2	176600836	176603836	NM_001177543_5_22.TTS-DW2000-5000	0	+
#chr2	176266768	176266848	NM_001177543_5_22.UTR5	0	+
#chr2	176270305	176270357	NM_001177543_5_22.UTR5	0	+
#chr2	176280240	176282881	NM_001177543_5_22.UTR3	0	+
#chr2	176598568	176598837	NM_001177543_5_22.UTR3	0	+
#chr11	70048970	70049197	NM_001177601_25.Intron.1	0	-
#chr11	70048905	70048970	NM_001177601_25.Coding_exon.1	0	-
#chr11	70049302	70049525	NM_001177601_25.Intron.2	0	-
#chr11	70049713	70050378	NM_001177601_25.Intron.3	0	-
#chr11	70050480	70051032	NM_001177601_25.Intron.4	0	-
#chr11	70051127	70051285	NM_001177601_25.Intron.5	0	-
#chr11	70049197	70049302	NM_001177601_25.Coding_exon.2	0	-
#chr11	70049525	70049713	NM_001177601_25.Coding_exon.3	0	-
#chr11	70050378	70050480	NM_001177601_25.Coding_exon.4	0	-
#chr11	70051032	70051127	NM_001177601_25.Coding_exon.5	0	-
#chr11	70051285	70051306	NM_001177601_25.Coding_exon.6	0	-
#chr11	70048705	70051416	NM_001177601_25	0610010K14Rik	-
#chr11	70048455	70048955	NM_001177601_25.TTS500	0	-
#chr11	70051165	70051665	NM_001177601_25.TSS500	0	-
#chr11	70048205	70049205	NM_001177601_25.TTS1000	0	-
#chr11	70050915	70051915	NM_001177601_25.TSS1000	0	-
#chr11	70047705	70049705	NM_001177601_25.TTS2000	0	-
#chr11	70050415	70052415	NM_001177601_25.TSS2000	0	-
#chr11	70046205	70051205	NM_001177601_25.TTS5000	0	-
#chr11	70048915	70053915	NM_001177601_25.TSS5000	0	-
#chr11	70048205	70048705	NM_001177601_25.TTS-DW0-500	0	-
#chr11	70051415	70051915	NM_001177601_25.TSS-UP0-500	0	-
#chr11	70047705	70048205	NM_001177601_25.TTS-DW500-1000	0	-
#chr11	70051915	70052415	NM_001177601_25.TSS-UP500-1000	0	-
#chr11	70046705	70047705	NM_001177601_25.TTS-DW1000-2000	0	-
#chr11	70052415	70053415	NM_001177601_25.TSS-UP1000-2000	0	-
#chr11	70043705	70046705	NM_001177601_25.TTS-DW2000-5000	0	-
#chr11	70053415	70056415	NM_001177601_25.TSS-UP2000-5000	0	-
#chr11	70048705	70048905	NM_001177601_25.UTR3	0	-
#chr11	70051306	70051416	NM_001177601_25.UTR5	0	-
#chr11	23475517	23476730	NM_027860_24.Intron.1	0	-
#chr11	23475480	23475517	NM_027860_24.Coding_exon.1	0	-
#chr11	23476822	23481595	NM_027860_24.Intron.2	0	-
#chr11	23481684	23482866	NM_027860_24.Intron.3	0	-
#chr11	23482940	23483097	NM_027860_24.Intron.4	0	-
#chr11	23483126	23484519	NM_027860_24.Intron.5	0	-
#chr11	23484611	23488634	NM_027860_24.Intron.6	0	-
#chr11	23488734	23489950	NM_027860_24.Intron.7	0	-
#chr11	23489983	23493425	NM_027860_24.Intron.8	0	-
#chr11	23493525	23495380	NM_027860_24.Intron.9	0	-
#chr11	23495526	23506613	NM_027860_24.Intron.10	0	-
#chr11	23506766	23509000	NM_027860_24.Intron.11	0	-
#chr11	23509124	23511746	NM_027860_24.Intron.12	0	-
#chr11	23511857	23511961	NM_027860_24.Intron.13	0	-
#chr11	23512048	23515095	NM_027860_24.Intron.14	0	-
#chr11	23515256	23517137	NM_027860_24.Intron.15	0	-
#chr11	23517196	23520228	NM_027860_24.Intron.16	0	-
#chr11	23520467	23522435	NM_027860_24.Intron.17	0	-
#chr11	23522529	23524935	NM_027860_24.Intron.18	0	-
#chr11	23525122	23528872	NM_027860_24.Intron.19	0	-
#chr11	23476730	23476822	NM_027860_24.Coding_exon.2	0	-
#chr11	23481595	23481684	NM_027860_24.Coding_exon.3	0	-
#chr11	23482866	23482940	NM_027860_24.Coding_exon.4	0	-
#chr11	23483097	23483126	NM_027860_24.Coding_exon.5	0	-
#chr11	23484519	23484611	NM_027860_24.Coding_exon.6	0	-
#chr11	23488634	23488734	NM_027860_24.Coding_exon.7	0	-
#chr11	23489950	23489983	NM_027860_24.Coding_exon.8	0	-
#chr11	23493425	23493525	NM_027860_24.Coding_exon.9	0	-
#chr11	23495380	23495526	NM_027860_24.Coding_exon.10	0	-
#chr11	23506613	23506766	NM_027860_24.Coding_exon.11	0	-
#chr11	23509000	23509124	NM_027860_24.Coding_exon.12	0	-
#chr11	23511746	23511857	NM_027860_24.Coding_exon.13	0	-
#chr11	23511961	23512048	NM_027860_24.Coding_exon.14	0	-
#chr11	23515095	23515256	NM_027860_24.Coding_exon.15	0	-
#chr11	23517137	23517196	NM_027860_24.Coding_exon.16	0	-
#chr11	23520228	23520467	NM_027860_24.Coding_exon.17	0	-
#chr11	23522435	23522529	NM_027860_24.Coding_exon.18	0	-
#chr11	23524935	23525122	NM_027860_24.Coding_exon.19	0	-
#chr11	23528872	23529022	NM_027860_24.Coding_exon.20	0	-
#chr11	23529031	23531204	NM_027860_24.Intron.20	0	-
#chr11	23531334	23533492	NM_027860_24.Intron.21	0	-
#chr11	23473775	23533631	NM_027860_24	0610010F05Rik	-
#chr11	23473525	23474025	NM_027860_24.TTS500	0	-
#chr11	23533380	23533880	NM_027860_24.TSS500	0	-
#chr11	23473275	23474275	NM_027860_24.TTS1000	0	-
#chr11	23533130	23534130	NM_027860_24.TSS1000	0	-
#chr11	23472775	23474775	NM_027860_24.TTS2000	0	-
#chr11	23532630	23534630	NM_027860_24.TSS2000	0	-
#chr11	23471275	23476275	NM_027860_24.TTS5000	0	-
#chr11	23531130	23536130	NM_027860_24.TSS5000	0	-
#chr11	23473275	23473775	NM_027860_24.TTS-DW0-500	0	-
#chr11	23533630	23534130	NM_027860_24.TSS-UP0-500	0	-
#chr11	23472775	23473275	NM_027860_24.TTS-DW500-1000	0	-
#chr11	23534130	23534630	NM_027860_24.TSS-UP500-1000	0	-
#chr11	23471775	23472775	NM_027860_24.TTS-DW1000-2000	0	-
#chr11	23534630	23535630	NM_027860_24.TSS-UP1000-2000	0	-
#chr11	23468775	23471775	NM_027860_24.TTS-DW2000-5000	0	-
#chr11	23535630	23538630	NM_027860_24.TSS-UP2000-5000	0	-
#chr11	23473775	23475480	NM_027860_24.UTR3	0	-
#chr11	23529022	23529031	NM_027860_24.UTR5	0	-
#chr11	23531204	23531334	NM_027860_24.UTR5	0	-
#chr11	23533492	23533631	NM_027860_24.UTR5	0	-
#chr5	65841676	65844424	NM_133697_61.Intron.1	0	-
#chr5	65841643	65841676	NM_133697_61.Coding_exon.1	0	-
#chr5	65844567	65851762	NM_133697_61.Intron.2	0	-
#chr5	65851811	65859631	NM_133697_61.Intron.3	0	-
#chr5	65844424	65844567	NM_133697_61.Coding_exon.2	0	-
#chr5	65851762	65851811	NM_133697_61.Coding_exon.3	0	-
#chr5	65859631	65859706	NM_133697_61.Coding_exon.4	0	-
#chr5	65859740	65883974	NM_133697_61.Intron.4	0	-
#chr5	65839993	65884074	NM_133697_61	1110003E01Rik	-
#chr5	65839743	65840243	NM_133697_61.TTS500	0	-
#chr5	65883823	65884323	NM_133697_61.TSS500	0	-
#chr5	65839493	65840493	NM_133697_61.TTS1000	0	-
#chr5	65883573	65884573	NM_133697_61.TSS1000	0	-
#chr5	65838993	65840993	NM_133697_61.TTS2000	0	-
#chr5	65883073	65885073	NM_133697_61.TSS2000	0	-
#chr5	65837493	65842493	NM_133697_61.TTS5000	0	-
#chr5	65881573	65886573	NM_133697_61.TSS5000	0	-
#chr5	65839493	65839993	NM_133697_61.TTS-DW0-500	0	-
#chr5	65884073	65884573	NM_133697_61.TSS-UP0-500	0	-
#chr5	65838993	65839493	NM_133697_61.TTS-DW500-1000	0	-
#chr5	65884573	65885073	NM_133697_61.TSS-UP500-1000	0	-
#chr5	65837993	65838993	NM_133697_61.TTS-DW1000-2000	0	-
#chr5	65885073	65886073	NM_133697_61.TSS-UP1000-2000	0	-
#chr5	65834993	65837993	NM_133697_61.TTS-DW2000-5000	0	-
#chr5	65886073	65889073	NM_133697_61.TSS-UP2000-5000	0	-
#chr5	65839993	65841643	NM_133697_61.UTR3	0	-
#chr5	65859706	65859740	NM_133697_61.UTR5	0	-
#chr5	65883974	65884074	NM_133697_61.UTR5	0	-
#'''
#'''test.bed.result---------------------------
#chr18	38409902	38410022	NM_024179_17.TcSS.dw200	0	+
#chr18	38410022	38410053	NM_024179_17.TcSS.dw200	0	+
#chr18	38410739	38410788	NM_024179_17.TcSS.dw200	0	+
#chr18	38409902	38410022	NM_024179_17.TcSS.dw400	0	+
#chr18	38410022	38410053	NM_024179_17.TcSS.dw400	0	+
#chr18	38410739	38410854	NM_024179_17.TcSS.dw400	0	+
#chr18	38410936	38411054	NM_024179_17.TcSS.dw400	0	+
#chr18	38412474	38412490	NM_024179_17.TcSS.dw400	0	+
#chr18	38420738	38420938	NM_024179_17.TsTS.both200	0	+
#chr18	38420938	38420939	NM_024179_17.TsTS.both200	0	+
#chr18	38421533	38421732	NM_024179_17.TsTS.both200	0	+
#chr18	38420717	38420938	NM_024179_17.TsTS.both400	0	+
#chr18	38419550	38419710	NM_024179_17.TsTS.both400	0	+
#chr18	38417997	38418016	NM_024179_17.TsTS.both400	0	+
#chr18	38420938	38420939	NM_024179_17.TsTS.both400	0	+
#chr18	38421533	38421932	NM_024179_17.TsTS.both400	0	+
#chr2	175163633	175163713	NM_001177543_18.TcSS.dw200	0	-
#chr2	175160126	175160178	NM_001177543_18.TcSS.dw200	0	-
#chr2	175160123	175160126	NM_001177543_18.TcSS.dw200	0	-
#chr2	175153034	175153099	NM_001177543_18.TcSS.dw200	0	-
#chr2	175160126	175160178	NM_001177543_18.TcSS.dw400	0	-
#chr2	175163633	175163713	NM_001177543_18.TcSS.dw400	0	-
#chr2	175160123	175160126	NM_001177543_18.TcSS.dw400	0	-
#chr2	175152972	175153099	NM_001177543_18.TcSS.dw400	0	-
#chr2	175152706	175152767	NM_001177543_18.TcSS.dw400	0	-
#chr2	175151472	175151549	NM_001177543_18.TcSS.dw400	0	-
#chr2	175150243	175150443	NM_001177543_18.TsTS.both200	0	-
#chr2	175150043	175150243	NM_001177543_18.TsTS.both200	0	-
#chr2	175150243	175150643	NM_001177543_18.TsTS.both400	0	-
#chr2	175149843	175150243	NM_001177543_18.TsTS.both400	0	-
#chr11	120209991	120210176	NR_038126_16.TcSS.dw200	0	+
#chr11	120212070	120212085	NR_038126_16.TcSS.dw200	0	+
#chr11	120209991	120210176	NR_038126_16.TcSS.dw400	0	+
#chr11	120212070	120212285	NR_038126_16.TcSS.dw400	0	+
#chr11	120212304	120212504	NR_038126_16.TsTS.both200	0	+
#chr11	120212104	120212504	NR_038126_16.TsTS.both400	0	+
#chr5	130695613	130695789	NM_001081394_4.TcSS.dw200	0	+
#chr5	130705318	130705337	NM_001081394_4.TcSS.dw200	0	+
#chr5	130705337	130705342	NM_001081394_4.TcSS.dw200	0	+
#chr5	130695613	130695789	NM_001081394_4.TcSS.dw400	0	+
#chr5	130705318	130705337	NM_001081394_4.TcSS.dw400	0	+
#chr5	130705337	130705496	NM_001081394_4.TcSS.dw400	0	+
#chr5	130707624	130707670	NM_001081394_4.TcSS.dw400	0	+
#chr5	130717147	130717168	NM_001081394_4.TsTS.both200	0	+
#chr5	130716217	130716361	NM_001081394_4.TsTS.both200	0	+
#chr5	130714265	130714300	NM_001081394_4.TsTS.both200	0	+
#chr5	130717168	130717368	NM_001081394_4.TsTS.both200	0	+
#chr5	130717147	130717168	NM_001081394_4.TsTS.both400	0	+
#chr5	130716217	130716361	NM_001081394_4.TsTS.both400	0	+
#chr5	130714116	130714300	NM_001081394_4.TsTS.both400	0	+
#chr5	130712772	130712823	NM_001081394_4.TsTS.both400	0	+
#chr5	130717168	130717568	NM_001081394_4.TsTS.both400	0	+
#chr7	52830496	52830546	NR_038165_1.TcSS.dw200	0	-
#chr7	52829997	52830147	NR_038165_1.TcSS.dw200	0	-
#chr7	52830496	52830546	NR_038165_1.TcSS.dw400	0	-
#chr7	52829977	52830147	NR_038165_1.TcSS.dw400	0	-
#chr7	52829782	52829892	NR_038165_1.TcSS.dw400	0	-
#chr7	52823679	52823749	NR_038165_1.TcSS.dw400	0	-
#chr7	52823164	52823364	NR_038165_1.TsTS.both200	0	-
#chr7	52823164	52823564	NR_038165_1.TsTS.both400	0	-
#chr7	52830496	52830546	NR_038166_2.TcSS.dw200	0	-
#chr7	52829997	52830147	NR_038166_2.TcSS.dw200	0	-
#chr7	52830496	52830546	NR_038166_2.TcSS.dw400	0	-
#chr7	52829977	52830147	NR_038166_2.TcSS.dw400	0	-
#chr7	52829782	52829892	NR_038166_2.TcSS.dw400	0	-
#chr7	52826492	52826562	NR_038166_2.TcSS.dw400	0	-
#chr7	52823164	52823364	NR_038166_2.TsTS.both200	0	-
#chr7	52823164	52823564	NR_038166_2.TsTS.both400	0	-
#chr5	31351012	31351045	NM_027855_3.TcSS.dw200	0	+
#chr5	31351045	31351129	NM_027855_3.TcSS.dw200	0	+
#chr5	31351834	31351917	NM_027855_3.TcSS.dw200	0	+
#chr5	31351012	31351045	NM_027855_3.TcSS.dw400	0	+
#chr5	31351045	31351129	NM_027855_3.TcSS.dw400	0	+
#chr5	31351834	31351953	NM_027855_3.TcSS.dw400	0	+
#chr5	31354569	31354641	NM_027855_3.TcSS.dw400	0	+
#chr5	31354834	31354906	NM_027855_3.TcSS.dw400	0	+
#chr5	31355135	31355155	NM_027855_3.TcSS.dw400	0	+
#chr5	31356635	31356740	NM_027855_3.TsTS.both200	0	+
#chr5	31356336	31356431	NM_027855_3.TsTS.both200	0	+
#chr5	31356740	31356940	NM_027855_3.TsTS.both200	0	+
#chr5	31356635	31356740	NM_027855_3.TsTS.both400	0	+
#chr5	31356333	31356431	NM_027855_3.TsTS.both400	0	+
#chr5	31355135	31355257	NM_027855_3.TsTS.both400	0	+
#chr5	31354834	31354906	NM_027855_3.TsTS.both400	0	+
#chr5	31354638	31354641	NM_027855_3.TsTS.both400	0	+
#chr5	31356740	31356996	NM_027855_3.TsTS.both400	0	+
#chr11	51502049	51502136	NM_025319_14.TcSS.dw200	0	-
#chr11	51499571	51499594	NM_025319_14.TcSS.dw200	0	-
#chr11	51499481	51499571	NM_025319_14.TcSS.dw200	0	-
#chr11	51499571	51499594	NM_025319_14.TcSS.dw400	0	-
#chr11	51502049	51502136	NM_025319_14.TcSS.dw400	0	-
#chr11	51499281	51499571	NM_025319_14.TcSS.dw400	0	-
#chr11	51499148	51499348	NM_025319_14.TsTS.both200	0	-
#chr11	51498948	51499148	NM_025319_14.TsTS.both200	0	-
#chr11	51499148	51499548	NM_025319_14.TsTS.both400	0	-
#chr11	51498886	51499148	NM_025319_14.TsTS.both400	0	-
#chr2	176266768	176266848	NM_001177543_5_22.TcSS.dw200	0	+
#chr2	176270305	176270357	NM_001177543_5_22.TcSS.dw200	0	+
#chr2	176270357	176270360	NM_001177543_5_22.TcSS.dw200	0	+
#chr2	176277384	176277449	NM_001177543_5_22.TcSS.dw200	0	+
#chr2	176266768	176266848	NM_001177543_5_22.TcSS.dw400	0	+
#chr2	176270305	176270357	NM_001177543_5_22.TcSS.dw400	0	+
#chr2	176270357	176270360	NM_001177543_5_22.TcSS.dw400	0	+
#chr2	176277384	176277511	NM_001177543_5_22.TcSS.dw400	0	+
#chr2	176277716	176277777	NM_001177543_5_22.TcSS.dw400	0	+
#chr2	176278934	176279011	NM_001177543_5_22.TcSS.dw400	0	+
#chr2	176280040	176280240	NM_001177543_5_22.TsTS.both200	0	+
#chr2	176280240	176280440	NM_001177543_5_22.TsTS.both200	0	+
#chr2	176279840	176280240	NM_001177543_5_22.TsTS.both400	0	+
#chr2	176280240	176280640	NM_001177543_5_22.TsTS.both400	0	+
#chr11	70051306	70051416	NM_001177601_25.TcSS.dw200	0	-
#chr11	70051285	70051306	NM_001177601_25.TcSS.dw200	0	-
#chr11	70051058	70051127	NM_001177601_25.TcSS.dw200	0	-
#chr11	70051306	70051416	NM_001177601_25.TcSS.dw400	0	-
#chr11	70051285	70051306	NM_001177601_25.TcSS.dw400	0	-
#chr11	70051032	70051127	NM_001177601_25.TcSS.dw400	0	-
#chr11	70050378	70050480	NM_001177601_25.TcSS.dw400	0	-
#chr11	70049641	70049713	NM_001177601_25.TcSS.dw400	0	-
#chr11	70048905	70048970	NM_001177601_25.TsTS.both200	0	-
#chr11	70049197	70049302	NM_001177601_25.TsTS.both200	0	-
#chr11	70049525	70049555	NM_001177601_25.TsTS.both200	0	-
#chr11	70048705	70048905	NM_001177601_25.TsTS.both200	0	-
#chr11	70048905	70048970	NM_001177601_25.TsTS.both400	0	-
#chr11	70049197	70049302	NM_001177601_25.TsTS.both400	0	-
#chr11	70049525	70049713	NM_001177601_25.TsTS.both400	0	-
#chr11	70050378	70050420	NM_001177601_25.TsTS.both400	0	-
#chr11	70048705	70048905	NM_001177601_25.TsTS.both400	0	-
#chr11	23533492	23533631	NM_027860_24.TcSS.dw200	0	-
#chr11	23531273	23531334	NM_027860_24.TcSS.dw200	0	-
#chr11	23529022	23529031	NM_027860_24.TcSS.dw400	0	-
#chr11	23531204	23531334	NM_027860_24.TcSS.dw400	0	-
#chr11	23533492	23533631	NM_027860_24.TcSS.dw400	0	-
#chr11	23528900	23529022	NM_027860_24.TcSS.dw400	0	-
#chr11	23475480	23475517	NM_027860_24.TsTS.both200	0	-
#chr11	23476730	23476822	NM_027860_24.TsTS.both200	0	-
#chr11	23481595	23481666	NM_027860_24.TsTS.both200	0	-
#chr11	23475280	23475480	NM_027860_24.TsTS.both200	0	-
#chr11	23475480	23475517	NM_027860_24.TsTS.both400	0	-
#chr11	23476730	23476822	NM_027860_24.TsTS.both400	0	-
#chr11	23481595	23481684	NM_027860_24.TsTS.both400	0	-
#chr11	23482866	23482940	NM_027860_24.TsTS.both400	0	-
#chr11	23483097	23483126	NM_027860_24.TsTS.both400	0	-
#chr11	23484519	23484598	NM_027860_24.TsTS.both400	0	-
#chr11	23475080	23475480	NM_027860_24.TsTS.both400	0	-
#chr5	65883974	65884074	NM_133697_61.TcSS.dw200	0	-
#chr5	65859706	65859740	NM_133697_61.TcSS.dw200	0	-
#chr5	65859640	65859706	NM_133697_61.TcSS.dw200	0	-
#chr5	65859706	65859740	NM_133697_61.TcSS.dw400	0	-
#chr5	65883974	65884074	NM_133697_61.TcSS.dw400	0	-
#chr5	65859631	65859706	NM_133697_61.TcSS.dw400	0	-
#chr5	65851762	65851811	NM_133697_61.TcSS.dw400	0	-
#chr5	65844425	65844567	NM_133697_61.TcSS.dw400	0	-
#chr5	65841643	65841676	NM_133697_61.TsTS.both200	0	-
#chr5	65844424	65844567	NM_133697_61.TsTS.both200	0	-
#chr5	65851762	65851786	NM_133697_61.TsTS.both200	0	-
#chr5	65841443	65841643	NM_133697_61.TsTS.both200	0	-
#chr5	65841643	65841676	NM_133697_61.TsTS.both400	0	-
#chr5	65844424	65844567	NM_133697_61.TsTS.both400	0	-
#chr5	65851762	65851811	NM_133697_61.TsTS.both400	0	-
#chr5	65859631	65859706	NM_133697_61.TsTS.both400	0	-
#chr5	65841243	65841643	NM_133697_61.TsTS.both400	0	-
#'''
