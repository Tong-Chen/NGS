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
    This is designed to count the number of predicted proteins in
    types like '3prime_partial', '5prime_partial', 'complete',
    'internal'.

    The output would be `file`.sta.xls and `file`.sta.xls.stackBars.pdf.
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
        metavar="FILEIN", help="Output of Transdecoder or \
parseTransdecoder.pep.py. Normally `Trinity.fasta.transdecoder.pep`.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_false", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    typeD = {}
    typeIDD = {}
    seqD  = {} 
    #type = ''
    for line in fh:
        if line[0] == '>':
            lineL = line.split()
            isoform = lineL[0].split('.')[0][1:]
            key = line
            #print lineL
            for item in lineL:
                #print item
                if item.startswith("type:"):
                    type = item.split(':')[1]
                    break
            if type not in typeD:
                typeD[type] = 0
                seqD[type] = 0
                typeIDD[type] = set()
            typeD[type] += 1
            typeIDD[type].add(isoform)
            #seqD[type][key] = 0
        else:
            seqD[type] += len(line.strip())-1
    #-------------------------------------------
    #print typeD
    output = file + '.sta.xls'
    outputID = file + '.id.xls'
    op_fh  = open(output, 'w')
    outputID_fh = open(outputID, 'w')
    
    print >>op_fh, "protein_type\tvalue\ttype"
    for type, cnt in typeD.items():
        print >>op_fh, "%s\t%s\t%s" % (type, cnt, "Protein count")

    for type, cnt in seqD.items():
        print >>op_fh, "%s\t%s\t%s" % (type, int(cnt/typeD[type]), "Protein mean length")

    op_fh.close()
    
    print >>outputID_fh, "ID\tType"
    for key, valueL in typeIDD.items():
        for isoform in valueL:
            print >>outputID_fh, "%s\t%s" % (isoform, key)

    cmd = "s-plot barPlot -f %s -m TRUE -a protein_type \
        -x 'Completeness of predicted proteins' -k free_y \
        -G protein_type -I value -P none -B type -O 1" % output
    os.system(cmd)
    
    
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


