#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright 2018, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to get promoter position from given bed file.

    The strand is considered.
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

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
        metavar="FILEIN", help="A bed file with at least six columns")
    parser.add_option("-u", "--upstream", dest="upstream",
        default='1000', help="Get upstream x bp as promoter. Default 1000.")
    parser.add_option("-d", "--downstream", dest="downstream",
        default='500', help="Get downstream x bp as promoter. Default 500.")
    parser.add_option("-n", "--addregiontoName", dest="addregiontoName",
        default=False, action="store_true", help="Add specified regions to names.")
    parser.add_option("-V", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    upstreamL = [int(i) for i in options.upstream.split(',')]
    downstreamL = [int(i) for i in options.downstream.split(',')]
    addregiontoName = options.addregiontoName
    
    if (len(upstreamL)>1):
        addregiontoName = True

    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    for line in fh:
        lineL = line.strip().split('\t')
        start = int(lineL[1])
        end   = int(lineL[2])
        name  = lineL[3]
        strand = lineL[5]
        tmpL = lineL[:]
        
        upL = []
        dwL = []
        if strand == '+':
            for up in upstreamL:
                p_up = start - up
                if(p_up) < 0:
                    p_up = 0
                upL.append(str(p_up))
            #---------------------
            for dw in downstreamL:
                p_dw = start + dw
                dwL.append(str(p_dw))
        elif strand == '-':
            for up in upstreamL:
                p_dw = end + up
                dwL.append(str(p_dw))
            for dw in downstreamL:
                p_up = end - dw
                if(p_up) < 0:
                    p_up = 0
                upL.append(str(p_up))
        else:
            print >>sys.stderr, "Unknown strand"
            sys.exit(1)

        for up, dw, up_sp, dw_sp in zip(upL, dwL, upstreamL, downstreamL):
            tmpL[1] = up
            tmpL[2] = dw
            if addregiontoName:
                tmpL[3] = name+':up'+str(up_sp)+'-dw'+str(dw_sp)
            print '\t'.join(tmpL)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
    if verbose:
        #print("--Successful %s" % strftime(timeformat, localtime()), file=sys.stderr)
        print >>sys.stderr, "--Successful %s" % strftime(timeformat, localtime())

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    with open('python.log', 'a') as fh:
        #print ("%s\n\tRun time : %s - %s " % \
        #(' '.join(sys.argv), startTime, endTime), file=fh)
        print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


