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
Input file format

PositionID      Offset  Sequence        Motif Name      Strand  MotifScore
224458  -303    GATTGCATTA  AARE(HLH)/mES-cMyc-ChIP-Seq/Homer   +       8.234019
224443  -102    TATTGCATCA  AARE(HLH)/mES-cMyc-ChIP-Seq/Homer   +       10.968466
225202  -423    AATTGCATCA  AARE(HLH)/mES-cMyc-ChIP-Seq/Homer   +       11.000554
225197  213     CATTGCATCA  AARE(HLH)/mES-cMyc-ChIP-Seq/Homer   +       10.816125
224953  -339    TATTGCATCA  AARE(HLH)/mES-cMyc-ChIP-Seq/Homer   -       9.4898

Output file format
Position    Motif1  Motif2  Motif3  Motif4  Motif5  Motif6  ......
224458  1   0   0   1   0   0
224443  0   2   6   0   0   0
225202  0   0   2   0   0   0
225197  3   0   0   1   0   0
224953  4   0   1   1   1   1
'''


import sys
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
def main():
    lensysargv = len(sys.argv)
    if lensysargv < 2:
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s filename \
column_for_rowname[1 based, default 1] \
column_for_colname[1 based, default 4] header[default 1]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    if lensysargv > 2:
        rowcol = int(sys.argv[2]) - 1
    else:
        rowcol = 0
    if lensysargv > 3:
        colcol = int(sys.argv[3]) - 1
    else:
        colcol = 3
    if lensysargv > 4:
        header = int(sys.argv[4])
    else:
        header = 1
    #-------------------------------------------------
    aDict = {}
    colname = set()
    for line in open(sys.argv[1]):
        if header:
            header -= 1
            continue
        #--------------------------------------------
        lineL = line.split()
        key = lineL[rowcol]
        value = lineL[colcol]
        colname.add(value)
        if key not in aDict:
            aDict[key] = {}
        if value not in aDict[key]:
            aDict[key][value] = 0
        aDict[key][value] += 1
    #----------END reading----------------------------
    colname = list(colname)
    colname.sort()
    print "peak\t%s" % "\t".join(colname)
    keyL = aDict.keys()
    keyL.sort()
    for key in keyL:
        tmpOutput = [key]
        for value in colname:
            if value in aDict[key]:
                tmpOutput.append(str(aDict[key][value]))
            else:
                tmpOutput.append('0')
        #-------------------------------------------
        print "\t".join(tmpOutput)
#------------------------------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


