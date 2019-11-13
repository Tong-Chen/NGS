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
    #Chr    Pos     Ref     Chain   Au_BS18_FM      Au_BS18_CD      Au_BS18_Diff

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
        metavar="FILEIN", help="")
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


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    group = options.group
    grpD = {}
    individualD = {}
    keepD = {}
    '''
    # Specify the sample difference status
    # H: high
    # L: low
    # N: no diff
    keepD = {"Au":['H', "L", "N"], "Nu":["H","L","N"]}
    '''
    verbose = options.verbose
    global debug
    debug = options.debug
    prefix = options.prefix
    #-----------------------------------
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
        keepD[grp] = set()
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
    diffSta = prefix + '.sta.xls'
    diffStaD = {}
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
            for i in range(4, len_lineL, 4):
                samp = '_'.join(lineL[i].split('_')[:2])
                assert samp in individualD, "Unknown sample %s " % samp
                grpIndexD[i] = individualD[samp]           
            header -= 1
            title = "%s\t%s" % (line, "_".join(grpL))
            continue
        #--------------------------------------
        keepD = {}
        for grp in grpL:
            keepD[grp] = []
        for i in range(4, len_lineL, 4):
            grp = grpIndexD[i]
            diff = float(lineL[i+2])
            p_value = float(lineL[i+3])
            type = 'N'
            if p_value <= 0.05:
                if diff <= -0.2:
                    type = 'L'
                elif diff >= 0.2:
                    type = 'H'
            #---------------------------------------
            keepD[grp].append(type)
        #---------------------------------------
        all = 1
        N_cnt = 0
        diffStaD_keyL = []
        for grp in grpL:
            typeL = keepD[grp]
            typeL.sort()
            diffStaD_keyL.append(''.join(typeL))
            if len(set(typeL)) > 1:
                all = 0
            else:
                if typeL[0] == "N":
                    N_cnt += 1
        descideDiff(diffStaD_keyL)
        diffStaD_key = tuple(diffStaD_keyL)
        diffStaD[diffStaD_key] = diffStaD.get(diffStaD_key, 0)+1
        if N_cnt == len(keepD):
            all = 0
        #------------------------------------------
        if debug:
            print >>sys.stderr, line
            for grp, typeL in keepD.items():
                print >>sys.stderr, grp, ''.join(typeL)
        if all:
            type = ''.join([keepD[grp][0] for grp in grpL])
            if type in outputD:
                fh = outputD[type]
            else:
                output = prefix + '.' + type + '.xls'
                fh = open(output, 'w')
                outputD[type] = fh
                print >>fh, title
            #---------------------------------------
            print >>fh, "%s\t%s" % (line, type)
    #-------------END reading file----------
    for fh in outputD.values():
        fh.close()
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    #--------diffSta output--------------------
    diffSta_fh = open(diffSta, 'w')
    print >>diffSta_fh, "%s\tCount" % '\t'.join(grpL)
    diffSta_keyL = diffStaD.keys()
    diffSta_keyL.sort()
    for diffSta_key in diffSta_keyL:
        print >>diffSta_fh, '%s\t%d' % ('\t'.join(diffSta_key), diffStaD[diffSta_key])
    diffSta_fh.close()
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


