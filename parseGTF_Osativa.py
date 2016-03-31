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
if False:
    print >>sys.stderr, "This program does not work under python 3, \
run in python 2.x."

from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"

'''
File format (The transcript IDs are assumed to be uniq especially for
neighboring genes. Lines containing CDS in the third column is not
necessary. Those lines will just be ignored. Only lines with 'exon',
'start_codon', 'stop_codon' are used.)


#gff-version 3
Chr1	MSU_osa1r7	gene	2903	10817	.	+	.	ID=LOC_Os01g01010;Name=LOC_Os01g01010;Note=TBC%20domain%20containing%20protein%2C%20expressed
Chr1	MSU_osa1r7	mRNA	2903	10817	.	+	.	ID=LOC_Os01g01010.1;Name=LOC_Os01g01010.1;Parent=LOC_Os01g01010
Chr1	MSU_osa1r7	exon	2903	3268	.	+	.	ID=LOC_Os01g01010.1:exon_1;Parent=LOC_Os01g01010.1
Chr1	MSU_osa1r7	exon	3354	3616	.	+	.	ID=LOC_Os01g01010.1:exon_2;Parent=LOC_Os01g01010.1
Chr1	MSU_osa1r7	exon	4357	4455	.	+	.	ID=LOC_Os01g01010.1:exon_3;Parent=LOC_Os01g01010.1
Chr1	MSU_osa1r7	exon	5457	5560	.	+	.	ID=LOC_Os01g01010.1:exon_4;Parent=LOC_Os01g01010.1
Chr1	MSU_osa1r7	exon	7136	7944	.	+	.	ID=LOC_Os01g01010.1:exon_5;Parent=LOC_Os01g01010.1
Chr1	MSU_osa1r7	exon	8028	8150	.	+	.	ID=LOC_Os01g01010.1:exon_6;Parent=LOC_Os01g01010.1
Chr1	MSU_osa1r7	exon	8232	8320	.	+	.	ID=LOC_Os01g01010.1:exon_7;Parent=LOC_Os01g01010.1

Tested for the following status (start_condon and stop_codon can span
two exons):
exon-ATT-exon-intron-exon-intron-exon-UAA-exon
exon-ATT-exon-UAA-exon
exon-AT-intron-T-exon-interon-exon-UA-intron-A-exon
exon-ex-ATT-on [means ATT at the last exon]

Attention:
    1.If only start_codon or stop_codon exists (not both ot all
    depleted), the start 3 letters or last 3 letters will be taken as
    the other codon.
    #the program may
    #output wrong results. You may want to complete the info about
    #start_codon or stop_codon. [positive strand with start codon and
    #negative strand with stopcodon, but without the other codon has
    #been solved.]

'''

def flankTSS_TES(TSS, TTS, end, regionL, strand, name_base, score, chr):
    #-------end means chromosome end-----------
    for seg in regionL:
        #--------------TSS---------------------------
        newstart = TSS - seg/2
        if newstart < 0: newstart = 0
        newend   = TSS + seg/2
        if newend > end: newend = end
        if strand == '+':
            name = name_base + '.TSS' + str(seg)
        elif strand == '-':
            name = name_base + '.TTS' + str(seg)
        else:
            name = name_base + 'unstranded' + str(seq)
        tmpLineL = [chr,str(newstart),str(newend),name, \
            score,strand]
        print "\t".join(tmpLineL)
        #------------TES---------------------------
        newstart = TTS - seg/2
        if newstart < 0: newstart = 0
        newend   = TTS + seg/2
        if newend > end: newend = end
        if strand == '+':
            name = name_base + '.TTS' + str(seg)
        elif strand == '-':
            name = name_base + '.TSS' + str(seg)
        tmpLineL = [chr,str(newstart),str(newend),name, \
            score,strand]
        print "\t".join(tmpLineL)
    #--------------------------------------
#------------------------------------------------------


def outTSS_TES(TSS, TTS, end, regionL, strand, name_base, score, chr):
    oldseg = 0
    for seg in regionL:
        #--------------TSS---------------------------
        newstart = TSS - seg
        if newstart < 0: newstart = 0
        newend = TSS - oldseg
        if newend > 0:
            if strand == '+':
                name = name_base + '.TSS-UP' + str(oldseg) + '-' + str(seg)
            elif strand == '-':
                name = name_base + '.TTS-DW' + str(oldseg) + '-' + str(seg)
            tmpLineL = [chr,str(newstart),str(newend),name, \
                score,strand]
            print "\t".join(tmpLineL)
        #------------TES---------------------------
        newstart = TTS + oldseg
        newend   = TTS + seg
        if newend > end: newend = end
        if newstart < end:
            if strand == '+':
                name = name_base + '.TTS-DW' + str(oldseg) + '-' + str(seg) 
            elif strand == '-':
                name = name_base + '.TSS-UP' + str(oldseg) + '-' + str(seg) 
            tmpLineL = [chr,str(newstart),str(newend),name, \
                score,strand]
            print "\t".join(tmpLineL)
        #--------update oldseg-------------------------
        oldseg = seg
    #------end for--------------------------------
#------------------------------------------------------

def get_gene_body(chr,keyL,tmpDict,name,score,strand):
    allpos_start = [ele[0] for ele in keyL]
    allpos_stop = [ele[1] for ele in keyL]
    allpos_start.sort()
    allpos_stop.sort()
    tmpLineL = [chr,str(allpos_start[0]),str(allpos_stop[-1]),name, score,strand]
    print '\t'.join(tmpLineL)
    return allpos_start[0], allpos_stop[-1]-1
#-------------------------------------------------------


def get_coding_exon_intron(chr,keyL,tmpDict,name,score,strand):
    oldEle = ''
    intron_num = 1
    exon_num = 1
    intron_start = -1
    begin_coding_exon = 0
    end_coding_exon   = 0
    oneMoreExonNeeded = 0
    ExonL = []
    tmpLineL_218 = ''
    first_met_stop_codon = 1
    first_met_start_codon = 1
    #print >>sys.stderr, set(tmpDict.values())
    if strand == '+' and \
            'start_codon' not in tmpDict.values() and \
            'stop_codon' in tmpDict.values():
        begin_coding_exon = 1
    elif strand == '-' and \
            'stop_codon' not in tmpDict.values() and \
            'start_codon' in tmpDict.values():
        begin_coding_exon = 1
    #print keyL
    #print >>sys.stderr, keyL
    #print >>sys.stderr, tmpDict
    for ele in keyL:
        if strand == '+':
            if tmpDict[ele] == 'exon':
                #--------intron----------------------
                if intron_start == -1:
                    intron_start = ele[1]
                else:
                    intron = name + '.Intron.'+str(intron_num)
                    intron_num += 1
                    tmpLineL = [chr,str(intron_start),str(ele[0]),intron,score,
                        strand]
                    intron_start = ele[1]
                    print '\t'.join(tmpLineL)
                #--------intron----------------------
                #---------output exon after start_codon meets---
                if begin_coding_exon and not end_coding_exon:
                    exon = name + '.Coding_exon.'+str(exon_num)
                    exon_num += 1
                    tmpLineL = [chr,str(ele[0]),str(ele[1]),exon,score,
                        strand]
                    ExonL.append(tmpLineL[:])
                    if tmpLineL_218:
                        print '\t'.join(tmpLineL_218)
                        tmpLineL_218 = ''
                    #print '\t'.join(tmpLineL)
                if oneMoreExonNeeded:
                    if exonEnd > ele[0]:
                        exon = name + '.Coding_exon.'+str(exon_num)
                        exon_num += 1
                        tmpLineL = [chr,str(ele[0]),str(exonEnd),exon,score,
                            strand]
                        print '\t'.join(tmpLineL)
                    oneMoreExonNeeded = 0
                #---------output exon after start_codon meets---
                oldEle = ele
                #-----------------------------------
            #------------------------------------------
            elif tmpDict[ele] == 'start_codon':
                if not begin_coding_exon:
                    #In case start_codon is the first one
                    if oldEle and oldEle[1] > ele[0] >= oldEle[0]:
                        exon = name + '.Coding_exon.'+str(exon_num)
                        exon_num += 1
                        tmpLineL_218 = [chr,str(ele[0]),str(oldEle[1]),exon,score,
                            strand]
                        #print '\t'.join(tmpLineL)
                    begin_coding_exon = 1
            elif tmpDict[ele] == 'stop_codon':
                #here = 0
                if not first_met_stop_codon:
                    exon = name + '.Coding_exon.'+str(exon_num)
                    exon_num += 1
                    tmpLineL = [chr,str(ele[0]),str(ele[1]),exon,score,
                        strand]
                    print '\t'.join(tmpLineL)
                elif first_met_stop_codon:
                    first_met_stop_codon = 0
                    for i in ExonL[:-1]:
                        #here = 1
                        print '\t'.join(i)
                    #if here:
                    #exon_num -= 1
                    #---------------If only one exon----------------
                    if tmpLineL_218 and ele[1] <= int(tmpLineL_218[2]):
                        tmpLineL_218[2] = str(ele[1])
                        print '\t'.join(tmpLineL_218)
                        tmpLineL_218 = ''
                    #---------------If only one exon----------------
                    elif tmpLineL_218 and ele[1] > int(tmpLineL_218[2]):
                        print '\t'.join(tmpLineL_218)
                        tmpLineL_218 = ''
                    elif oldEle[0] < ele[1] <= oldEle[1]:
                        exon_num -= 1
                        exon = name + '.Coding_exon.'+str(exon_num)
                        exon_num += 1
                        tmpLineL = [chr,str(oldEle[0]),str(ele[1]),exon,score,
                            strand]
                        print '\t'.join(tmpLineL)
                    if ele[1] > oldEle[1]:
                        #print oldEle, ExonL
                        if ExonL:
                            assert int(ExonL[-1][2]) == oldEle[1]
                            print '\t'.join(ExonL[-1])
                        oneMoreExonNeeded = 1
                        exonEnd = ele[1]
                    #---this can never happen----
                    #if ele[1] > oldEle[1]:
                    #    oneMoreExonNeeded = 1
                    #    exonEnd = ele[1]
                    end_coding_exon = 1
            #--------END +----------------------------------------
        #-------------END + strand----------------------
        elif strand == '-':
            #print >>sys.stderr, exon_num
            if tmpDict[ele] == 'exon':
                #---------output intron-------------------------
                if intron_start == -1:
                    intron_start = ele[1]
                else:
                    intron = name + '.Intron.'+str(intron_num)
                    intron_num += 1
                    tmpLineL = [chr,str(intron_start),str(ele[0]),intron,score,
                        strand]
                    intron_start = ele[1]
                    print '\t'.join(tmpLineL)
                #---------output intron-------------------------
                #---------output exon after start_codon meets---
                if begin_coding_exon and not end_coding_exon:
                    exon = name + '.Coding_exon.'+str(exon_num)
                    exon_num += 1
                    tmpLineL = [chr,str(ele[0]),str(ele[1]),exon,score,
                        strand]
                    ExonL.append(tmpLineL[:])
                    #print >>sys.stderr, ExonL
                    #print >>sys.stderr, ExonL
                    #---output exon defined by stop_codon
                    if tmpLineL_218:
                        print '\t'.join(tmpLineL_218)
                        tmpLineL_218 = ''
                    new_exon = 1
                    #print '\t'.join(tmpLineL)
                if oneMoreExonNeeded:
                    if exonEnd > ele[0]:
                        exon = name + '.Coding_exon.'+str(exon_num)
                        exon_num += 1
                        tmpLineL = [chr,str(ele[0]),str(exonEnd),exon,score,
                            strand]
                        print '\t'.join(tmpLineL)
                    oneMoreExonNeeded = 0
                #---------output exon after start_codon meets---
                oldEle = ele
                #print >>sys.stderr, oldEle
            #----------------------END if--------------------
            elif tmpDict[ele] == 'stop_codon':
                if not begin_coding_exon:
                    #In case stop_codon is the first one
                    if oldEle and oldEle[1] > ele[0] >= oldEle[0]:
                        exon = name + '.Coding_exon.'+str(exon_num)
                        exon_num += 1
                        tmpLineL_218 = [chr,str(ele[0]),str(oldEle[1]),exon,score,
                            strand]
                        #print >>sys.stderr, '\t'.join(tmpLineL_218)
                    begin_coding_exon = 1
                    new_exon = 0
                    #print >>sys.stderr, oldEle
                    #print >>sys.stderr, ele
            elif tmpDict[ele] == 'start_codon':
                #here = 0
                if not first_met_start_codon:
                    #if tmpLineL_218 and ele[1] > int(tmpLineL_218[2]):
                        #----two exons,  but start_codon is the begin
                        #of last exon
                    #    print '\t'.join(tmpLineL_218)
                    #    tmpLineL_218 = ''
                    exon = name + '.Coding_exon.'+str(exon_num)
                    exon_num += 1
                    tmpLineL = [chr,str(ele[0]),str(ele[1]),exon,score,
                        strand]
                    print '\t'.join(tmpLineL)
                    #oneMoreExonNeeded = 0
                elif first_met_start_codon:
                    #print ele[1],oldEle[1],tmpLineL_218
                    first_met_start_codon = 0
                    #print >>sys.stderr, "here"
                    #print >>sys.stderr, ExonL
                    for i in ExonL[:-1]:
                        #here = 1
                        print '\t'.join(i)
                    #if here:
                    #exon_num -= 1
                    #---------------If only one exon----------------
                    if tmpLineL_218 and ele[1] <= int(tmpLineL_218[2]):
                        tmpLineL_218[2] = str(ele[1])
                        print '\t'.join(tmpLineL_218)
                        tmpLineL_218 = ''
                    #---------------If only one exon----------------
                    #This will hit----------------------
                    elif tmpLineL_218 and ele[1] > int(tmpLineL_218[2]):
                        #----two exons,  but start_codon is the begin
                        #of last exon
                        print '\t'.join(tmpLineL_218)
                        tmpLineL_218 = ''
                    elif oldEle[0] < ele[1] <= oldEle[1]:
                        exon_num -= 1
                        exon = name + '.Coding_exon.'+str(exon_num)
                        exon_num += 1
                        tmpLineL = [chr,str(oldEle[0]),str(ele[1]),exon,score,
                            strand]
                        print '\t'.join(tmpLineL)
                    #This will hit----------------------
                    if ele[1] > oldEle[1]:
                        #--Output only when new exons met--
                        #print oldEle, ExonL
                        if ExonL:
                            assert int(ExonL[-1][2]) == oldEle[1]
                            print '\t'.join(ExonL[-1])
                        #--Output only when new exons met--
                        oneMoreExonNeeded = 1
                        exonEnd = ele[1]
                    end_coding_exon = 1
        #--------END -----------------------------------------
    #------------END for ele------------------------------
    #----for NR--------------
    if first_met_stop_codon and first_met_start_codon:
        #---lose one side of codon--
        if begin_coding_exon:
            #print >>sys.stderr, 'Here'
            if tmpLineL_218:
                print '\t'.join(tmpLineL_218)
                tmpLineL_218 = ''
            for i in ExonL[:]:
                print '\t'.join(i)
        else:
            for ele in keyL:
                exon = name + '.Coding_exon.'+str(exon_num)
                exon_num += 1
                tmpLineL = [chr,str(ele[0]),str(ele[1]),exon,score,
                    strand]
                print '\t'.join(tmpLineL)
    #----for NR--------------
#-------------------end get_coding_exon_intron_genebody--------
'''
tmpDict = {strand:'+',  (0, 100):'exon', (80, 83):'start_codon', \
        (1000, 1003):'stop_codon', (700, 1500):'exon'}
'''

def parse(tmpDict, key, name_label,cs_dict,regionL):
    strand = tmpDict.pop('strand')
    chr    = tmpDict.pop('chr')
    score, name, unused_315 = key.split('#')
    name += '_'+name_label
    keyL = tmpDict.keys()
    #keyL.sort(key=lambda x:x[0])
    keyL.sort()
    get_coding_exon_intron(chr,keyL,tmpDict,name,'0',strand)
    start, end = get_gene_body(chr,keyL,tmpDict,name,score,strand)
    if cs_dict:
        chrom_end = cs_dict[chr]
        flankTSS_TES(start, end, chrom_end, regionL, strand, name, '0', chr)
        outTSS_TES(start, end, chrom_end, regionL, strand, name, '0', chr)
    if strand == '+':
        #------find UTR5--------------
        tmpL = []
        needOneMoreExon = 0
        all_UTR3 = 0
        #---stop_codons span two exons---
        stop_codon_unfinished = 0
        start_codon_unfinished = 0
        stop_codon_cnt = 0
        for key in keyL:
            if tmpDict[key] == "exon":
                tmpL.append(key)
                if all_UTR3 and not needOneMoreExon:
                    tmpLineL = [chr,str(key[0]),str(key[1]),UTR3,'0',
                        strand]
                    print '\t'.join(tmpLineL)
                    continue
            elif tmpDict[key] == "start_codon":
                #if start_codon is the first one, no UTR5
                if not tmpL:
                    continue
                #if start codon spans multile exons, ignore its end
                #position
                if start_codon_unfinished:
                    continue
                UTR5 = name + '.UTR5'
                #------for UTRs span multiple exons----
                #intron_start = -1
                if len(tmpL) > 1:
                    for ele in tmpL[:-1]:
                        tmpLineL = [chr,str(ele[0]),str(ele[1]),UTR5,'0',
                            strand]
                        print '\t'.join(tmpLineL)
                        #------------intron-------------------------
                        #if intron_start == -1:
                        #    intron_start = ele[1]
                        #else:
                        #    intron = name + '.Intron.'+str(intron_num)
                        #    intron_num += 1
                        #    tmpLineL = [chr,str(intron_start),str(ele[0]),intron,'0',
                        #       strand]
                        #    intron_start = ele[1]
                        #    print '\t'.join(tmpLineL)
                        #-----------intron-------------------------
                #-----for UTRs overlapped with exon---------------------- 
                #--Only need the start site of start_condon-----
                #----------------intron----------------------------
                #if intron_start == -1:
                #    intron_start = tmpL[-1][1]
                #else:
                #    intron = name + '.Intron.'+str(intron_num)
                #    intron_num += 1
                #    tmpLineL = [chr,intron_start,tmpL[-1][0],intron,'0',
                #        strand]
                #    intron_start = tmpL[-1][1]
                #    print '\t'.join(tmpLineL)
                #----------------intron----------------------------
                #This assert can not be right if the start_codon begins
                #at the next exon.
                #assert tmpL[-1][0] <= key[0] <= tmpL[-1][1]
                #assert tmpL[-1][0] <= key[0] <= tmpL[-1][1]
                if tmpL[-1][1] >= key[0] > tmpL[-1][0]:
                    tmpLineL = [chr,str(tmpL[-1][0]),str(key[0]),UTR5,'0',
                        strand]
                    print '\t'.join(tmpLineL)
                    #-------------coding exon---------------------------
                    #if key[0] < tmpL[-1][1]:
                    #    cod_exon = name + '.Coding_exon.' + str(exon_num)
                    #    exon_num += 1
                    #    tmpLineL = [chr,key[0],tmpL[-1][1],cod_exon,score,
                    #        strand]
                    #    print '\t'.join(tmpLineL)
                    #-------------coding exon---------------------------
                elif key[0] > tmpL[-1][1]:
                    tmpLineL = [chr,str(tmpL[-1][0]),str(tmpL[-1][1]),UTR5,'0',
                        strand]
                    print '\t'.join(tmpLineL)
                if key[1] - key[0] < 3:
                    start_codon_unfinished = 1
            #-------------------end UTR5------------------------
            elif tmpDict[key] == "stop_codon":
                #--this is used to record the end of stop_codon
                if stop_codon_unfinished:
                    end_stop_codon = key[1]
                    stop_codon_cnt += key[1]-key[0]
                    #this continue is usd to read in another exon
                #--this is used to record the end of stop_codon
                else:
                    all_UTR3 = 1
                    UTR3 = name + '.UTR3'
                    lastExon = tmpL[-1]
                    #This assert could not be allowed when stop codon
                    #begins at one exon.
                    #assert lastExon[0] <= key[0] <= lastExon[1], name
                    ##if key[1] == lastExon[1],  UTR must start at the
                    #next exon 
                    if key[1] < lastExon[1]:
                        tmpLineL = [chr,str(key[1]),str(lastExon[1]),UTR3,'0',
                            strand]
                        print '\t'.join(tmpLineL)
                    elif key[1] == lastExon[1] and key[1]-key[0]<3:
                        stop_codon_unfinished = 1
                        stop_codon_cnt += key[1]-key[0]
                        needOneMoreExon = 1
                        end_stop_codon  = key[1]
                        tmpL = []
                    elif key[1] > lastExon[1] and key[0] < lastExon[1]:
                        print >>sys.stderr, "This should never happen"
                        sys.exit(1)
                        needOneMoreExon = 1
                        end_stop_codon  = key[1]
                        tmpL = []
                        stop_codon_cnt = 3
                    if key[0] >= lastExon[1]:
                        needOneMoreExon = 1
                        stop_codon_cnt = key[1]-key[0]
                        end_stop_codon = key[1]
                        tmpL = []
                    #------------------------------
                #---------------------------------
            #---------------------------------------
            if needOneMoreExon and tmpL and stop_codon_cnt == 3:
                assert (len(tmpL) == 1) 
                if end_stop_codon >= tmpL[0][0]:
                    tmpLineL = [chr,str(end_stop_codon),str(tmpL[0][1]),UTR3,'0',
                        strand]
                    print '\t'.join(tmpLineL)
                else:
                    tmpLineL = [chr,str(tmpL[0][0]),str(tmpL[0][1]),UTR3,'0',
                        strand]
                    print '\t'.join(tmpLineL)
                #--------------------------------------
                needOneMoreExon = 0           
            #--------End <needOneMoreExon>-------------
        #-------------end iterating---------------
        #-------------find UTR5 and UTR3----------------
    elif strand == '-':
        #TSS = keyL[-1][1] - 1
        #TTS = keyL[0][0]
        #flankTSS_TES(TSS, TTS, chrom_end, regionL, strand, name, '0', chr)
        #outTSS_TES(TSS, TTS, chrom_end, regionL, strand, name, '0', chr)
        #------find UTR5--------------
        tmpL = []
        needOneMoreExon = 0
        all_UTR5 = 0
        start_codon_unfinished = 0
        stop_codon_unfinished = 0
        start_codon_cnt = 0
        #print keyL
        for key in keyL:
            if tmpDict[key] == "exon":
                tmpL.append(key)
                if all_UTR5 and not needOneMoreExon:
                    tmpLineL = [chr,str(key[0]),str(key[1]),UTR5,'0',
                        strand]
                    print '\t'.join(tmpLineL)
                    continue
            elif tmpDict[key] == "stop_codon":
                #if stop_codon is the first one, no UTR3
                if not tmpL:
                    continue
                if stop_codon_unfinished:
                    continue
                UTR3 = name + '.UTR3'
                #------for UTRs span multiple exons
                if len(tmpL) > 1:
                    for ele in tmpL[:-1]:
                        tmpLineL = [chr,str(ele[0]),str(ele[1]),UTR3,'0',
                            strand]
                        print '\t'.join(tmpLineL)
                #-----for UTRs overlapped with exon---------------------- 
                #--Only need the start site of start_condon-----
                #This assert can not be right if the stop_codon begins
                #at the next exon.
                #assert tmpL[-1][0] <= key[0] <= tmpL[-1][1]
                if tmpL[-1][1] >= key[0] > tmpL[-1][0]:
                    tmpLineL = [chr,str(tmpL[-1][0]),str(key[0]),UTR3,'0',
                        strand]
                    print '\t'.join(tmpLineL)
                elif key[0] > tmpL[-1][1]:
                    tmpLineL = [chr,str(tmpL[-1][0]),str(tmpL[-1][1]),UTR3,'0',
                        strand]
                    print '\t'.join(tmpLineL)
                if key[1] - key[0] < 3:
                    stop_codon_unfinished = 1
            #-------------------end UTR3------------------------
            elif tmpDict[key] == "start_codon":
                #--this is used to record the end of stop_codon
                if start_codon_unfinished:
                    begin_start_codon = key[1]
                    start_codon_cnt += key[1]-key[0]
                    #this continue is usd to read in another exon
                    #continue
                #--this is used to record the end of stop_codon
                else:
                    all_UTR5 = 1
                    UTR5 = name + '.UTR5'
                    #print tmpL
                    lastExon = tmpL[-1]
                    #This assert could not be allowed when start codon
                    #begins at one exon.
                    #assert lastExon[0] <= key[0] <= lastExon[1], name
                    #This assert is used to make sure if codon spans
                    #multiple exons,  it will be separated into several
                    #segments.
                    #assert key[1] <= lastExon[1]
                    if key[1] < lastExon[1]:
                        tmpLineL = [chr,str(key[1]),str(lastExon[1]),UTR5,'0',
                            strand]
                        print '\t'.join(tmpLineL)
                    elif key[1] == lastExon[1] and key[1]-key[0]<3:
                        start_codon_unfinished = 1
                        start_codon_cnt += key[1]-key[0]
                        needOneMoreExon = 1
                        begin_start_codon  = key[1]
                        tmpL = []
                    elif key[1] > lastExon[1] and key[0] < lastExon[1]:
                        print >>sys.stderr, key[1],lastExon[1],name
                        print >>sys.stderr, "This should never happen"
                        sys.exit(1)
                        needOneMoreExon = 1
                        begin_start_codon = key[1]
                        tmpL = []
                        start_codon_cnt = key[1]-key[0]
                    if key[0] >= lastExon[1]:
                        needOneMoreExon = 1
                        start_codon_cnt = key[1]-key[0]
                        begin_start_codon = key[1]
                        tmpL = []
                    #---------------------------------
                #------------------------------------
            #----to deal with stop_codon at the end of exon--------------------
            if needOneMoreExon and tmpL and start_codon_cnt == 3:
                if begin_start_codon >= tmpL[0][0]:
                    tmpLineL = [chr,str(begin_start_codon),str(tmpL[0][1]),UTR5,'0',
                        strand]
                    print '\t'.join(tmpLineL)
                else:
                    tmpLineL = [chr,str(tmpL[0][0]),str(tmpL[0][1]),UTR5,'0',
                        strand]
                    print '\t'.join(tmpLineL)
                #--------------------------------------
                needOneMoreExon = 0           
            #--------End <needOneMoreExon>-------------
        #-------------end iterating---------------
        #-------------find UTR5 and UTR3----------------
#---------------------------------------

def main():
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, "***Print results to screen.***"
        print >>sys.stderr, "This parses GTF file to extract UTR5, \
UTR3, coding exon, intron, TSS, TTS, upstream and downstream regions."
        print >>sys.stderr, 'Using python %s filename[GTF, sorted \
first by gene, then column 4 and 5(using program sortGTF.py); - means \
sys.stdin] chromsome_size [if given, upstream and downstream regions \
will be output; Otherwise only outputting inner regions]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    regionL = [500,1000,2000,5000]
    accepted_regionL = ['exon', 'stop_codon', 'start_codon']
    #------get chromosome size-----------------
    chrom_size = ''
    cs_dict = {}
    if len(sys.argv) > 2:
        chrom_size = sys.argv[2]
        for line in open(chrom_size):
            lineL = line.split()
            cs_dict[lineL[0]] = int(lineL[1])
    #------get chromosome size-----------------
    file = sys.argv[1]
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #-------------------------------------------
    name_label = 0
    tmpDict = {}
    oldname = ""
    for line in fh:
        lineL = line.strip().split("\t")
        region_name = lineL[2]
        # GTF coordinate starts with 1. Here transfer to bed format.
        if region_name in accepted_regionL:
            start_end = (int(lineL[3])-1,int(lineL[4])) 
            #print >>sys.stderr, line
            #print >>sys.stderr, lineL
            attributeL = lineL[8].split(';')
            #print >>sys.stderr, attributeL
            #-----------------Get gene_id---------------------
            transcript_id = attributeL[0].replace("ID=","")
            gene_id = ''
            for attribute in attributeL[1:]:
                if attribute.find("Parent=") != -1:
                    gene_id = attribute.replace("Parent=", "")
                    break
            if not gene_id:
                for attribute in attributeL[1:]:
                    if attribute.find("Name=") != -1:
                        gene_id = attribute.replace("Name=", "")
                        break
            #-----------------Get gene_id---------------------
            key = '#'.join([gene_id, transcript_id, lineL[0]])
            if oldname != "" and oldname != key:
                name_label += 1
                parse(tmpDict[oldname],oldname,str(name_label),cs_dict,regionL)
                tmpDict = {}
            #-------------------------------------------------------
            oldname = key
            if key not in tmpDict:
                tmpDict[key] = {}
                tmpDict[key]['strand'] = lineL[6]
                tmpDict[key]['chr'] = lineL[0]
            if start_end not in tmpDict[key]:
                tmpDict[key][start_end] = region_name
            else:
                print >>sys.stderr, "Duplicate start_end (%d, %d) \
                    for %s" % (start_end[0], start_end[1], key)
                sys.exit(1)
            #--------------------------------------------------
        #--------------end filtering lines-----------------------
    #------------------end reading----------------------------
    #---deal with the last element in tmpDict-----------
    if tmpDict:
        name_label += 1
        parse(tmpDict[oldname],oldname,str(name_label),cs_dict,regionL)
        tmpDict = {}
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
#------------------------------------------------------
    #-------------END reading file----------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


