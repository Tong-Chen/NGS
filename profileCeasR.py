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
#import re
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"

#main = re.compilr('main="(.+?)"')

def readCeasR(file, adict):
    lineL = []
    for line in open(file):
        if line.startswith('pie(x=x,labels=c('):
            if line.find('Genome') != -1:
                key = 'pie_Genome_' + file
                if key not in adict:
                    adict[key] = []
            elif line.find('ChIP') != -1:
                key = 'pie_ChIP_' + file
                if key not in adict:
                    adict[key] = []
        elif line.startswith('legend("top",legend=c("Promoter'):
            start = line.find('=c("') + 4
            end = line.find('"),col')
            lineL = line[start:end].replace(' %','').split('","')
            adict[key] = lineL[:]
            key = ''
        elif line.find('Average Profile near TSS') != -1:
            lineL.append(line)
            key = 'TSS_' + file
            adict[key] = lineL[-3:]
            key = ''
        elif line.find('Average Profile near TTS') != -1:
            key = 'TTS_' + file
            adict[key] = lineL[-2:]
            key = ''
        elif line.find('Average Gene Profile') != -1:
            key = 'Gene_' + file
            adict[key] = lineL[-2:]
            key = ''
        else:
            lineL.append(line.strip())
#--------------------------------
def outputPie(adict, keyL, prefix):
    file = prefix + '.pie'
    fh = open(file, 'w')
    regL = ''
    for item in keyL:
        if item.startswith('pie_Genome'):
            if not regL:
                regL = [pos.split(':')[0].strip() for pos in adict[item]]
                posL = [pos.split(':')[1].strip() for pos in adict[item]]
                print >>fh, "#Sample\t%s" % "\t".join(regL)
                print >>fh, "Genome\t%s" % '\t'.join(posL)
        elif item.startswith('pie_ChIP_'):
            #regL = [pos.split(':')[0].strip() for pos in adict[item]]
            posL = [pos.split(':')[1].strip() for pos in adict[item]]
            #print "Sample\t%s" % "\t".join(regL)
            print >>fh, "%s\t%s" % (item.split('_',2)[2] ,'\t'.join(posL))
        #-------------------------------------------------------------
    #--------------------------------------------------
    fh.close()
#---------------------------------------------------------------------
def outputProfile(adict, keyL, type, prefix):
    file = prefix + '.' + type
    fh = open(file, 'w')
    posL = []
    newKeyL = []
    first = 1
    for item in keyL:
        if item.startswith(type):
            if first:
                posL.append(adict[item][0][5:-1].split(','))
                newKeyL.append('POS')
                first = 0
            #---------------------------------------
            posL.append(adict[item][1][5:-1].split(','))
            newKeyL.append(item.split('_',1)[1])
    #----------------------------------------------------------
    print >>fh, '\t'.join(newKeyL)
    if posL:
        lenPos = len(posL[0])
        for i in range(0, lenPos):
            print >>fh, '\t'.join([posLItem[i] for posLItem in posL])
    #---------------------------------------------
    fh.close()
#--------------------------------------------------
def main():
    if len(sys.argv) < 3:
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s prefix ceaseR1 ceaseR2...' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    prefix = sys.argv[1]
    adict = {}
    for file in sys.argv[2:]:
        readCeasR(file, adict)
    #--------------------------------
    keyL = adict.keys()
    keyL.sort(reverse=True)
    #keyL = sys.argv[2:]
    #print adict
    outputPie(adict, keyL, prefix)
    outputProfile(adict, keyL, 'TSS', prefix)
    outputProfile(adict, keyL, 'TTS', prefix)
    outputProfile(adict, keyL, 'Gene', prefix)
    
#------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    #startTime = localtime()
    #startTime = '-'.join([str(x) for x in startTime[:3]]) \
    #    + ' ' + ':'.join([str(x) for x in startTime[3:6]])
    main()
    endTime = strftime(timeformat, localtime())
    #endTime = localtime()
    #endTime = '-'.join([str(x) for x in endTime[:3]]) \
    #    + ' ' + ':'.join([str(x) for x in endTime[3:6]])
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


