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
#UTR and Coding_exon information is not necessary unless you want
#to separate junctions in coding region or uncoding regions.  
chr5	31351012	31351045	NM_027855_3.UTR5	0   +
chr5	31356740	31356996	NM_027855_3.UTR3	0   +
chr5	130695613	130695789	NM_001081394_4.UTR5	0   +
chr5	130705318	130705337	NM_001081394_41081394_4.UTR5    0	+
chr5	130717168	130719635	NM_001081394_4.UTR3	0   +
chr3	152090270	152091796	NM_198416_30606.Coding_exon.1	0   +
chr3	152109762	152109901	NM_198416_30606.Coding_exon.2   0	+
chr3	152111742	152111890	NM_198416_30606.Coding_exon.3   0	+
chr3	152111976	152112092	NM_198416_30606.Coding_exon.4   0	+
chr3	1521126802112608	152112680	NM_198416_30606.Coding_exon.5   0	+
'''


import sys
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
def main():
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, "This program is designed to get \
junction regions from a bed file."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename ' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    typeL = ['UTR5','UTR3','Coding_exon']
    lenL = [100,50,30]
    aDict = {}
    #aDict = {gene:{UTR5:[(pos1,pos2),()],UTR3:[],Coding_exon:[]}}
    for line in open(sys.argv[1]):
        lineL = line.split()
        keyL  = lineL[3].split('.')
        key   = keyL[0]
        type  = keyL[1]
        if key not in aDict:
            aDict[key] = {}
        if type not in aDict[key]:
            aDict[key][type] = []
        aDict[key][type].append(lineL)
    #-------------------------------------------------
    for gene,itemD in aDict.items():
        juncD = {}
        #juncD = {100:[[(),()]],50:[[(),()]]}
        tmpL = []
        UTR5_n = 0
        CE_n   = 0 #Coding_exon
        UTR3_n = 0
        if 'UTR5' in itemD:
            UTR5L = itemD['UTR5']
            tmpL.extend(UTR5L)
            UTR5_n = len(UTR5L)
        if 'UTR3' in itemD:
            UTR3L = itemD['UTR3']
            tmpL.extend(UTR3L)
            UTR3_n = len(UTR3L)
        if 'Coding_exon' in itemD:
            CEL = itemD['Coding_exon']
            tmpL.extend(CEL)
            CE_n = len(CEL)
        #-------------------------------------------------
        tmpL.sort(key=lambda x:(int(x[1]),int(x[2])))
        #print >>sys.stderr,tmpL
        outputL = tmpL[0]
        start1 = int(tmpL[0][1])
        end1   = int(tmpL[0][2])
        len1   = end1 - start1
        chr    = tmpL[0][0]
        strand = tmpL[0][-1]
        for i in tmpL[1:]:
            start2 = int(i[1])
            end2   = int(i[2])
            len2   = end2 - start2
            #print >>sys.stderr,end1,start2
            #print >>sys.stderr,juncD
            if end1 < start2:
                for len_t in lenL:
                    if len_t not in juncD:
                        juncD[len_t] = []
                    if len1 >=2*len_t and len2 >=2*len_t:
                        juncD[len_t].append([(end1-len_t,end1),\
                            (start2,start2+len_t)])           
                    else:
                        juncD[len_t].append([(-1,-1),(-1,-1)])
                #---------------------------------
                start1 = start2
                end1   = end2
                len1   = len2
            elif end1 == start2:
                for len_t in lenL:
                    if len_t not in juncD:
                        juncD[len_t] = []
                    juncD[len_t].append([(-2,-2),(-2,-2)])
                #----------------------------------
                end1 = end2
                len1 = end1 - start1
            else:
                #print >>sys.stderr, "This shoudl never happend"
                sys.exit(1)
            #------------------------------------------
        #----get all juncs for one isoform--------------------
        for key,itemL in juncD.items():
            if strand == '-':
                itemL.reverse()   
            #--------------------------------
            cnt = 0
            #print >>sys.stderr,cnt
            for i in range(UTR5_n,1,-1):
                junc = itemL[cnt]
                cnt += 1
                for i124 in range(2):
                    outputL[1] = str(junc[i124][0])
                    outputL[2] = str(junc[i124][1])
                    outputL[3] = '@'.join([gene,'UTR5','J',str(key),str(cnt)])
                    print "\t".join(outputL)
                #print >>sys.stderr,cnt
            if UTR5_n > 0:
                junc = itemL[cnt]
                cnt += 1
                for i124 in range(2):
                    outputL[1] = str(junc[i124][0])
                    outputL[2] = str(junc[i124][1])
                    outputL[3] = '@'.join([gene,'UTR5:CE','J',str(key),str(cnt)])
                    print "\t".join(outputL)
                #print >>sys.stderr,cnt
            for i in range(CE_n,1,-1):
                junc = itemL[cnt]
                cnt += 1
                for i124 in range(2):
                    outputL[1] = str(junc[i124][0])
                    outputL[2] = str(junc[i124][1])
                    outputL[3] = '@'.join([gene,'CE','J',str(key),str(cnt)])
                    print "\t".join(outputL)
                #print >>sys.stderr,cnt
            if CE_n > 0 and UTR3_n > 0:
                junc = itemL[cnt]
                cnt += 1
                for i124 in range(2):
                    outputL[1] = str(junc[i124][0])
                    outputL[2] = str(junc[i124][1])
                    outputL[3] = '@'.join([gene,'CE:UTR3','J',str(key),str(cnt)])
                    print "\t".join(outputL)
                #print >>sys.stderr,cnt
            for i in range(UTR3_n,1,-1):
                junc = itemL[cnt]
                cnt += 1
                for i124 in range(2):
                    outputL[1] = str(junc[i124][0])
                    outputL[2] = str(junc[i124][1])
                    outputL[3] = '@'.join([gene,'UTR3','J',str(key),str(cnt)])
                    print "\t".join(outputL)
                #print >>sys.stderr,cnt
            #----------Finishe output one type of len-----------------
        #---------------Finish output one refseq---------------------
    #-------------------Finish all isofomrs----

#------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


