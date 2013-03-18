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
    if lensysargv != 2:
        print >>sys.stderr, "Get hcp, icp and lowcp promoters \
according to paper doi:10.1038/nature06008]"
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s fastafile' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    print "Gene\toeCpG\tfracGC\tType"
    for line in open(sys.argv[1]):
        if line[0] == '>':
            key = line[1:-1]
        else:
            interval = 500
            slide = 1
            len_line = len(line) - 1
            tmpL = []
            for i in range(0,len_line-interval+1,slide):
                assert(i+interval <= len_line)
                segment = line[i:i+interval]
                CpG = segment.count('CG')
                c = segment.count('C')
                g = segment.count('G')
                if (c*g != 0):
                    oeCpG = (CpG * 1.0 / (c * g)) * interval
                else:
                    oeCpG = 0

                fracGC = (c+g) * 1.0 / interval
                #print "%s\t%f\t%f\t%f\t%f\t%f" % (segment, CpG,
                #        c,g,oeCpG, fracGC)
                tmpL.append((oeCpG,fracGC))
            #-------------------------------------------
            tmpL.sort(key=lambda x: x[0], reverse=True)
            #print tmpL
            type = "ICP"
            if tmpL[0][0] <= 0.4:
                type = "LCP"
            else:
                for item in tmpL:
                    if item[0] < 0.6:
                        break
                    if item[0] >=0.6 and item[1] >= 0.55:
                        type = "HCP"
                        break
            #----------------------------------------------
            print "%s\t%.1f\t%.1f\t%s" % \
                (key, tmpL[0][0], tmpL[0][1], type)
        #------------------------------------------
    #------------------------------------------
#----------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


