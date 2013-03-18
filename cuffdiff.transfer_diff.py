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
def main():
    lensysargv = len(sys.argv)
    if lensysargv < 3:
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename \
output_prefix least_rpkm[1] least_length[300]' % sys.argv[0]
        sys.exit(0)
    #-----parameter processing-------------------------
    file = sys.argv[1]
    prefix = sys.argv[2]
    if lensysargv > 3:
        least_rpkm = float(sys.argv[3])
    else:
        least_rpkm = 1
    if lensysargv > 4:
        least_length = int(sys.argv[4])
    else:
        least_length = 300
    #------------------------------------------------------
    aDict = {}
    header = 1
    for line in open(file):
        if header:
            header -= 1
            header_four = line.split()[:4]
            #header_four = line.split()[:1]
            continue
        #--------------------------------------------
        lineL = line.split()
        key = '\t'.join(lineL[:4])
        samp1,samp2,status,v1,v2,sig = \
            lineL[4],lineL[5],lineL[6],lineL[7],lineL[8],lineL[-1]
        #key = lineL[0]
        #samp1,samp2,v1,v2,sig = \
        #    lineL[1],lineL[2],lineL[4],lineL[5],lineL[-1]
        if key not in aDict:
            aDict[key] = dict()
        if samp1 not in aDict[key]:
            aDict[key][samp1] = v1
        else:
            if status == 'OK':
                assert abs(float(v1)-float(aDict[key][samp1]))<1 ,key
        if samp2 not in aDict[key]:
            aDict[key][samp2] = v2
        else:
            if status == 'OK':
                assert abs(float(v2)-float(aDict[key][samp2]))<1 ,key
            #assert v2 == aDict[key][samp2], key
        sig_key = samp1+'___'+samp2
        if sig_key not in aDict[key]:
            aDict[key][sig_key] = sig
        else:
            print "Duplicate %s" % key
            sys.exit(1)
    #-----------END reading---------------------------------
    #------------output header---------------
    file = prefix + '.table'
    fh = open(file, 'w')
    file_expr = prefix + '.expr_' + str(least_rpkm)
    fh_expr = open(file_expr, 'w')
    file_expr_len = prefix + '.expr_' + str(least_rpkm) \
        + '.len_' + str(least_length)
    fh_expr_len = open(file_expr_len, 'w')
    headerL = aDict[key].keys()
    headerL.sort(key=lambda x: sum([ord(y) for y in list(x)]))
    sampL = [x for x in headerL if x.find('___') == -1]
    combL = [x for x in headerL if x.find('___') != -1]
    #----------output header------------------------------
    print >>fh, "%s\t%s" % ('\t'.join(header_four), '\t'.join(headerL))
    print >>fh_expr, "%s\t%s" % ('\t'.join(header_four), '\t'.join(headerL))
    print >>fh_expr_len, "%s\t%s" % ('\t'.join(header_four), '\t'.join(headerL))
    keyL = aDict.keys()
    keyL.sort()
    for key in keyL:
        valueL = [aDict[key][innerK] for innerK in headerL]
        print >>fh, "%s\t%s" % (key, '\t'.join(valueL))
        #------output genes which have rpkm at least 1 in 
        #------at least ones sample --------------------
        output = 0
        for samp in sampL:
            if float(aDict[key][samp]) >= least_rpkm:
                output = 1
                break
        #------output transcripts with reasonable length--------------
        if output:
            print >>fh_expr, "%s\t%s" % (key, '\t'.join(valueL))
            start_end = key.split('\t')[-1].split(':')[-1].split('-')
            if int(start_end[1]) - int(start_end[0]) + 1 \
                >= least_length:
                print >>fh_expr_len, "%s\t%s" % (key, '\t'.join(valueL))
        #--------------------------------------------------
    fh.close()
    fh_expr.close()
    fh_expr_len.close()
    #----------------------------------------------------
    #------output expressed ones------------------------    
#----------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


