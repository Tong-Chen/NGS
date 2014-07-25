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
desc = '''
Functional description:
    This is designed to filter the result of
    <relativePositonOfMultipleSites.py>.
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="The output of \
relativePositonOfMultipleSites.py.")
    parser.add_option("-s", "--sep", dest="sep",
        default=':', help="The separtor between peak_name and \
other things in the third column. Default ':'.")
    parser.add_option("-c", "--compare-col-begin", dest="colB",
        default=13, help="The start column containing compare \
information, simply the column containing @. A zero-started \
number needed. Default 13, means the forth column.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    sep = options.sep
    colB = int(options.colB)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    peakD = {}
    '''
    peakD = {peak: {
                    key_target_region1: 
                        { 
                            posType1: [posDist, line]
                            posType2: [posDist, line]
                        }
                    key_target_region2: 
                        { 
                            posType1: [posDist, line]
                            posType2: [posDist, line]
                        }
            }
    '''
    targetTypeL = set()
    posTypeL = set()
    for line in fh:
        if line.find("No motif") == -1 and \
            line.find("Untargeted") == -1:
            lineL = line.split()
            peak = lineL[3].split(sep)[0]
            if peak not in peakD:
                peakD[peak] = {}
            for relativePos in lineL[colB:]:
                pos, key = relativePos.split('@')
                targetTypeL.add(key)
                posType, posDist = pos.split("_")
                posTypeL.add(posType)
                if key not in peakD[peak]:
                    peakD[peak][key] = {}
                #-------------------------------------
                if posType not in peakD[peak][key]:
                    peakD[peak][key][posType] = [int(posDist), line]
                else:
                    if posType in ["UP", "DW"]:
                        if posDist < peakD[peak][key][posType][0]:
                            peakD[peak][key][posType] = [int(posDist), line]
                        #--------------------------------
                    elif posType in ["OUP", "ODW"]:
                        if posDist > peakD[peak][key][posType][0]:
                            peakD[peak][key][posType] = [int(posDist), line]
                    elif posType == "IN":
                        #No operation
                        pass
                #-----------------------------------------------------
            #-----------------------------------------
        #---------------------------------------------
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    #--------------prepare output -----------------
    '''
    peakD = {peak: {
                    key_target_region1: 
                        { 
                            posType1: [posDist, line], 
                            posType2: [posDist, line]
                        }
                    key_target_region2: 
                        { 
                            posType1: [posDist, line], 
                            posType2: [posDist, line]
                        }
            }
    '''
    print "peak\t%s" % '\t'.join(targetTypeL)
    for peak, targetD in peakD.items():
        tmpL = []
        for targetType in targetTypeL:
            relativePosD = targetD[targetType]
            if "IN" in relativePosD:
                tmp = "IN"
            elif "OUP" in relativePosD or "ODW" in relativePosD:
                if "OUP" not in relativePosD:
                    ODW_dist = relativePosD["ODW"][0]
                    tmp = "ODW_"+str(ODW_dist)
                elif "ODW" not in relativePosD:
                    OUP_dist = relativePosD["OUP"][0]
                    tmp = "OUP_"+str(OUP_dist)
                else:
                    OUP_dist = relativePosD["OUP"][0]
                    ODW_dist = relativePosD["ODW"][0]
                    tmp = "OUP_"+str(OUP_dist) \
                        if OUP_dist <=ODW_dist else \
                         "ODW_"+str(ODW_dist)
                #-------------------------------------------
            elif "UP" in relativePosD or "DW" in relativePosD:
                if "UP" not in relativePosD:
                    DW_dist = relativePosD["DW"][0]
                    tmp = "DW_"+str(DW_dist)
                elif "DW" not in relativePosD:
                    UP_dist = relativePosD["UP"][0]
                    tmp = "UP_"+str(UP_dist)
                else:
                    UP_dist = relativePosD["UP"][0]
                    DW_dist = relativePosD["DW"][0]
                    tmp = "UP_"+str(UP_dist) \
                        if UP_dist <=DW_dist else \
                         "DW_"+str(DW_dist)
                #------------------------------------
            #----------------------------
            tmpL.append(tmp)
        #---------------------------------------
        print "%s\t%s" % (peak, '\t'.join(tmpL))
    #---------------------------------------------------------
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



