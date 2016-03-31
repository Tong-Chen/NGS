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

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from bs4 import BeautifulSoup

#reload(sys)
#sys.setdefaultencoding('utf8')

#from multiprocessing.dummy import Pool as ThreadPool

debug = 0

def fprint(content):
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
        metavar="FILEIN", help="The output of \
Trinotate/Trinotate.sqlite report.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def process(item, label, tmpL, header=0):
    blastL = ['sprot_Top_BLASTX_hit', 'TrEMBL_Top_BLASTX_hit',
        'sprot_Top_BLASTP_hit', 'TrEMBL_Top_BLASTP_hit']
    #'Pfam', 'SignalP', 'TmHMM', 'eggnog',
    eggnog = ['eggnog']
    pfam = ['Pfam']
    #go = ['gene_ontology_blast', 'gene_ontology_pfam']
    depleted = ['RNAMMER', 'prot_coords', 'transcript', 'peptide']
    directReturn = ['prot_id', 'gene_ontology_blast', 'gene_ontology_pfam']
    if label in blastL:
        if header:
            tmpL.append(item+'_id')
            tmpL.append(item+'_desp')
        else:
            if item == '.':
                tmpL.append('.')
                tmpL.append('.')
            else:
                itemL = item.split('^')
                tmpL.append(itemL[0])
                tmpL.append(itemL[5].replace('Full=',\
                    '').split('{')[0].\
                    replace('RecName: ', "").replace("SubName: ", ""))
    elif label in directReturn:
        tmpL.append(item)
    #elif label in go:
    #    if header:
    #        tmpL.append(item)
    #    else:
    #        tmpL.append(item)
    elif label in pfam:
        if header:
            #tmpL.append(item+'_id')
            #tmpL.append(item+'_name')
            #tmpL.append(item+'_desp')
            tmpL.append(item)
        else:
            if item == '.':
                tmpL.append('.')
            else:
                itemL = item.split('`')
                valueTl = []
                for item in itemL:
                    valueTl.append('^'.join(item.split('^')[:3]))
                tmpL.append('`'.join(valueTl))
                #tmpL.append(itemL[0])
                #tmpL.append(itemL[1])
                #tmpL.append(itemL[2])
    elif label in depleted:
        pass
#---------------------------------

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
    header = 1
    for line in fh:
        tmpL = []
        if header:
            headerL = line.strip().split('\t')
            tmpL.extend(headerL[:2])
            len_header = len(headerL)
            header -= 1
            for i in range(2, len_header):
                item  = headerL[i]
                label = headerL[i]
                process(item, label, tmpL, 1)
            print '\t'.join(tmpL)
            continue
        lineL = line.strip().split('\t')
        tmpL.extend(lineL[:2])
        otuputL = lineL[:2]
        for i in range(2, len_header):
            item  = lineL[i]
            label = headerL[i]
            process(item, label, tmpL)
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


