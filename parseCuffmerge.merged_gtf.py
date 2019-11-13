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
if False:
    print "This program does not work under python 3, \
run in python 2.x."

import sys
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
def main():
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, "This file is used to count number of \
various types from merged.gtf, the output of cuffmerge. \
Print the result to two files."
        print >>sys.stderr, 'Using python %s filename[- means \
sys.stdin]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    classCode = { \
        'u':'Unknown_intergenic', \
        'i':'Intron_derived', \
        'j':'Novel_isoforms', \
        'x':'Exon_antisense', \
        'o':'Exonic_overlap', \
        'c':'Contained', \
        'p':'Run-on_frag', \
        '=':'Known', \
        'e':'Pre-mRNA', \
        'r':'Repeats', \
        's':'Intron_antisense', \
        '.':'Multiple_classification'
        }

    definedKeyL = ['gene_id', 'transcript_id',\
        'gene_name','oLd', 'nearest_ref', 'class_code']
    file = sys.argv[1]
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #-------------------------
    geneDict = {}
    trDict   = {}
    oldTr = ''
    header = 0
    withRef = ['=', 'j', 's', 'x', 'o', 'c', 'i', 'p', 'e', '.']
    newISO = ['u', 'r' ]
    tmpL = set()
    for line in fh:
        if header:
            header -= 1
            continue
        #here is your reading
        #print line
        if oldTr == '' or line.find(oldTr) == -1: 
            lineL = line.strip().strip(';').split('\t')
            focus = lineL[-1]
            focusDict = dict([[i.split(" ")[0], \
                i.split(" ")[1].replace('"', '')] \
                for i in focus.split("; ")])
            #print >>sys.stderr, focus
            #print >>sys.stderr, focusDict
            #tmpL = tmpL.union(set(focusDict.keys()))
            #print focusDict
            class_code = focusDict['class_code']
            if 'contained_in' in focusDict:
                class_code = 'c'
            if class_code in [".", "-"]:
                if 'gene_name' not in focusDict:
                    class_code = 'u'
            gn = focusDict['gene_id']
            tr = focusDict['transcript_id']
            oldTr = 'transcript_id "'+tr+'"'
            #oldTr = '"'+tr+'"'
            if gn not in geneDict:
                geneDict[gn] = []
                if class_code in withRef:
                    geneDict[gn].append([focusDict['gene_name'], gn, class_code])
                    #----------------------------------------
                elif class_code in newISO:
                    geneDict[gn].append([\
                        '.'.join(focusDict['oId'].split('.')[:-1]),gn, class_code])
                else:
                    print >>sys.stderr, line
                    print >>sys.stderr, focusDict
                    print >>sys.stderr, "UNpected class_code \
<%s> for <%s>" % (class_code, tr)
                    sys.exit(1)
            else:
                if class_code in withRef:
                    gene_name = focusDict['gene_name']
                    new = 1
                    for itemL in geneDict[gn]:
                        if gene_name == itemL[0]:
                            new = 0
                            break
                    #-----------------------------
                    if new:
                        geneDict[gn].append([gene_name,gn,class_code])
                    #appear = tmp89L[0].rfind(gene_name)
                    #if appear ==-1 or \
                    #(len(tmp89L[0])>appear+len(gene_name) and \
                    #tmp89L[0][appear+len(gene_name)] != '&'):
                    #    tmp89L = geneDict[gn].split('__')
                    #    tmp89L[0] += '&'+focusDict['gene_name']
                    #    geneDict[gn] = '__'.join(tmp89L)
                elif class_code in newISO:
                    oId = '.'.join(focusDict['oId'].split('.')[:-1])
                    new = 1
                    for itemL in geneDict[gn]:
                        if oId == itemL[0]:
                            new = 0
                            break
                    #-----------------------------
                    if new:
                        geneDict[gn].append([oId,gn,class_code])
                    #appear = tmp89L[0].rfind(oId)
                    #if appear==-1 or \
                    #(len(tmp89L[0])>appear+len(oId) and \
                    #tmp89L[0][appear+len(oId)] != '&'):
                    #    tmp89L = geneDict[gn].split('__')
                    #    tmp89L[0] += '&'+oId
                    #    geneDict[gn] = '__'.join(tmp89L)
                else:
                    print >>sys.stderr, line
                    print >>sys.stderr, focusDict
                    print >>sys.stderr, "UNpected class_situation \
<%s> for <%s>" % (class_code, tr)
                    sys.exit(1)
            #--------------------------------------------------
            if tr not in trDict:
                if class_code in withRef:
                    trDict[tr] = '__'.join([focusDict['gene_name'], \
                        focusDict['nearest_ref'], tr, class_code])
                elif class_code in newISO:
                    trDict[tr] = '__'.join([tr, class_code])
                else:
                    print >>sys.stderr, line
                    print >>sys.stderr, focusDict
                    print >>sys.stderr, "UNpected class_code <%s> \
for <%s>" % (class_code, tr)
                    sys.exit(1)
            else:
                print >>sys.stderr, "Should not be here %s" % tr
                sys.exit(1)
            #--aEND = nd j ------------------------------------------
        #-------------------------------------------
    #-------------END reading file----------
    #print tmpL
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    output = open(file+'.gname', 'w')
    geneKeyL = geneDict.keys()
    geneKeyL.sort()
    for gene in geneKeyL:
        for itemL in geneDict[gene]:
            #print itemL
            print >>output, "%s\t%s" % (gene, '__'.join(itemL))
    output.close()

    output = open(file+'.trname', 'w')
    trKeyL = trDict.keys()
    trKeyL.sort()
    for tr in trKeyL:
        print >>output, "%s\t%s" % (tr, trDict[tr])
    output.close()

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


