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

def main():
    if len(sys.argv) != 5:
        print >>sys.stderr, '''
  This script is used to parse fimo output. First it will get the
  occurance percentage of each motif. Second it will get the position
  of the motif relative to summit and real position in genome.
    Print the result to two files,  output_prefix.sta &
    output_prefix.bed.

'''
        print >>sys.stderr, 'Using python %s fimo.txt peak.bed \
summits.bed output_prefix' % sys.argv[0]
        sys.exit(0)
    #----------------------------------------------
    fimo = sys.argv[1]
    fimoD = {}
    peak = sys.argv[2]
    peakD = dict([(line.split()[3], line.split()[0:3]) \
            for line in open(peak)])
    summit = sys.argv[3]
    summitD = dict([(line.split()[3], line.split()[1]) \
            for line in open(summit)])
    output_prefix = sys.argv[4]
    #-----------------------------------------------------
    head = 1
    for line in open(fimo):
        if head:
            head -= 1
            continue
        #--------------------------
        motif, peak, start, stop = line.split('\t', 4)[:4]
            #-------------position
        start = int(start)
        stop = int(stop)
        if start > stop:
            start, stop = stop, start
        realStart = int(peakD[peak][1]) + start - 1
        realStop  = int(peakD[peak][1]) + stop
        length_p = int(peakD[peak][2]) - int(peakD[peak][1])
            #--------relative position-----
        realMedian = (stop+start)/2.0 + int(peakD[peak][1]) - 1
        relativeDist = realMedian - int(summitD[peak]) 
        percentageDist = relativeDist / length_p
        posStr = ','.join([str(realStart), str(realStop),
            str(relativeDist)]) 
        if motif not in fimoD:
            fimoD[motif] = {}
        if peak not in fimoD[motif]:
            fimoD[motif][peak] = [peakD[peak][0], str(realStart),
                    str(realStop), relativeDist, percentageDist, posStr]
        else:
            if relativeDist < 0:
                relativeDist_tmp = (-1) * relativeDist
            else:
                relativeDist_tmp = relativeDist
            #-------------------------------------
            old_relativeDist = fimoD[motif][peak][3]
            if old_relativeDist < 0:
                old_relativeDist_tmp = (-1) * old_relativeDist
            else:
                old_relativeDist_tmp = old_relativeDist
            #--------------------------------------
            if relativeDist_tmp < old_relativeDist_tmp:
                fimoD[motif][peak][1] = str(realStart)
                fimoD[motif][peak][2] = str(realStop)
                fimoD[motif][peak][3] = relativeDist
                fimoD[motif][peak][4] = percentageDist
            #----------------------------
            fimoD[motif][peak].append(posStr)
        #---------------------END one line-------------------------
    #-------------------------END-----------------
    #------------output-------------------
    '''
    fimo = {motif:{peak:[chr, start, end, relativeDist, posStr...]}}
    '''
    fh1 = open(output_prefix+'.sta', 'w')
    fh2 = open(output_prefix+'.bed', 'w')
    for key, valueD in fimoD.items():
        print >>fh1, '%s\t%d' % (key, len(valueD))       
        for peak, bedL in valueD.items():
            print >>fh2, "%s\t%s\t%s\t%s\t%.1f\t%s\t%.2f\t%s" % \
                (bedL[0], bedL[1], bedL[2], key, bedL[3], peak,
                       bedL[4], ';'.join(bedL[5:]))
    fh1.close()
    fh2.close()
#----------------------------------------------------------
if __name__ == '__main__':
    main()

