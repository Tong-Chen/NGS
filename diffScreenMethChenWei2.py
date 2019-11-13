#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to screen out diffMethC.

file:
 #Chr Pos Ref Chain Au_BS18_F Au_BS18_M Au_BS18_C Au_BS18_D Au_E101_F Au_E101_M

group:
    Au_BS18 Au
    Au_E101 Au

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from scipy import stats
from math import log

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

def fprint(content):
    """ 
    This is a Google style docs.

    Args:
        param1(str): this is the first param
        param2(int, optional): this is a second param
            
    Returns:
        bool: This is a description of what is returned
            
    Raises:
        KeyError: raises an exception))
    """
    print json_dumps(content,indent=1)

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
        metavar="FILEIN", help="Merged FMCD file")
    parser.add_option("-g", "--groupFile", dest="group",
        help="A group file with first column as individual and second column as group. First line will be treated as header line and will be skipped.")
    parser.add_option("-o", "--output-prefix", dest="prefix",
        help="Prefix for output files.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def descideDiff(diffStaD_keyL):
    '''
    This is designed to decide which type of result should be output.

    diffStaD_keyL 
        = ['HHH', 'HH'] -- output
    or
        = ['HHN', 'LL'] -- output
    or
        = ['HHH', 'NN'] -- output
    '''
    pass
#-------------------------------------------
def meanDiff(valueL1, valueL2, add_value=0.05):
    meanV1 = sum(valueL1) / len(valueL1)
    meanV2 = sum(valueL2) / len(valueL2)
    return meanV1-meanV2
    #return log((meanV1+add_value) / (meanV2+add_value))/log(2)
#----------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    group = options.group
    grpD = {}
    individualD = {}
    verbose = options.verbose
    global debug
    debug = options.debug
    prefix = options.prefix
    #-----------------------------------
    ttest = stats.ttest_ind
    group_header = 1
    for line in open(group):
        if group_header:
            group_header -= 1
            continue
        individual, grp = line.split()
        individualD[individual] = grp
        if grp not in grpD:
            grpD[grp] = [individual]
        else:
            grpD[grp].append(individual)
    #=---------------------------
    #print >>sys.stderr, individualD
    grpL = grpD.keys()
    grpL.sort()
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    outputD = {}
    """
    diffStaD = {('NNN', 'LLL'):5, ('HHH', 'LHN'):6}
    """
    #diffSta = prefix + '.sta.xls'
    #diffStaD = {}
    header = 1
    grpIndexD = {}
    '''
    #Map position to group
    grpIndexD = {4: "Au", 8:"Au", 12:"Nu"}
    '''
    for line in fh:
        line = line.strip()
        lineL = line.split('\t')
        len_lineL = len(lineL)
        #keyL = lineL[:4]
        if header:
            for i in range(4, len_lineL):
                samp, fmcd = lineL[i].split('@')
                assert samp in individualD, "Unknown sample %s " % samp
                grpIndexD[i] = individualD[samp]           
            header -= 1
            title = "%s\t%s\t%s" % (line, "-".join(grpL), "p_value")
            continue
        #--------------------------------------
        diffD = {}
        '''
        diffD = {'grp1':[1, 2, 3], 'grp2':[4, 5, 6]}
        '''
        for i in range(4, len_lineL):
            grp = grpIndexD[i]
            if grp not in diffD:
                diffD[grp] = [float(lineL[i])]
            else:
                diffD[grp].append(float(lineL[i]))
        #---------------------------------------
        all =0
        grp0 = grpL[0]
        grp1 = grpL[1]
        value0 = diffD[grp0]
        value1 = diffD[grp1]
        valueDiff = meanDiff(value0, value1)
        if abs(valueDiff) < 0.001:
            p_value = 1
        else:
            p_value = ttest(value0, value1)[1]
        if debug:
            print >>sys.stderr, line
            print >>sys.stderr, valueDiff, p_value
        if p_value <= 0.05:
            if valueDiff < 0:
                type = grp1+'_up'
                all = 1
            elif valueDiff > 0:
                type = grp0+'_up'
                all = 1
        #------------------------------------------
        #if debug:
        #    print >>sys.stderr, line
        #    for grp, typeL in keepD.items():
        #        print >>sys.stderr, grp, ''.join(typeL)
        if all:
            if type in outputD:
                fh = outputD[type]
            else:
                output = prefix + '.' + type + '.xls'
                fh = open(output, 'w')
                outputD[type] = fh
                print >>fh, title
            #---------------------------------------
            print >>fh, "%s\t%.2f\t%.4f" % (line, valueDiff, p_value)
    #-------------END reading file----------
    for fh in outputD.values():
        fh.close()
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    #--------diffSta output--------------------
    #diffSta_fh = open(diffSta, 'w')
    #print >>diffSta_fh, "%s\tCount" % '\t'.join(grpL)
    #diffSta_keyL = diffStaD.keys()
    #diffSta_keyL.sort()
    #for diffSta_key in diffSta_keyL:
    #    print >>diffSta_fh, '%s\t%d' % ('\t'.join(diffSta_key), diffStaD[diffSta_key])
    #diffSta_fh.close()
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
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
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


