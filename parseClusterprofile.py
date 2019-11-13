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
    This is designed to parse the output of `clusterProfile.sh`, KEGG
    enrichment file. 

Normally <input file format>

ID Description GeneRatio BgRatio pvalue p.adjust qvalue geneID Count
hsa05200 Pathways in cancer 113/1176 397/7082 7e-10 2e-7 1e-7 7046/5727/4313/2776/5568/999/2736 113

The program will search for <geneID> or any other string given to <-s> 
in the first line to find the column containing gene ids.

One can also supply a number like `8` to <-c> to tell the program the
genes ids are in the eighth column.


The gene ids are `entrez id`, a map file containing both gene names
and ids should be supplied.

<ID_map file format> <tab> delimier for columns
Symbol	EntrezID
A1BG	1
A1BG-AS1	503538
A1CF	29974
A2M	2
A2M-AS1	144571
A2ML1	144568
A2MP1	3
A3GALT2	127550
A4GALT	53947
A4GNT	51146
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
        metavar="FILEIN", help="Normally the output of \
`clusterProfile.sh` as format specified above.")
    parser.add_option("-s", "--string-for-geneID-column",
        dest="col_name", default='geneID', 
        help="The name of entrez gene ID column. Default <geneID>.")
    parser.add_option("-c", "--column-index-for-geneID", dest="col_index",
        type="int", help="A number like `8` to specify gene ids \
are contained in the 8th column. This parameter has large priority \
than <-s>. No default value for this parameter.")
    parser.add_option("-m", "--map-id", dest="idMap",
        help="ID map file as specified above, with the first column \
as gene names and the second column as entrez gene ids.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readIdMap(idMap):
    idMapD = {}
    for line in open(idMap):
        gene, id = line.strip().split('\t')
        assert id not in idMapD, "Ambigious id %s" % id
        idMapD[id] = gene
    return idMapD
#-------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    col_name  = options.col_name
    col_index = options.col_index
    idMap     = options.idMap
    idMapD    = readIdMap(idMap)
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    file_out = file + '.update'
    file_out_fh = open(file_out, 'w')
    header = 1
    for line in fh:
        lineL = line.strip().split('\t')
        lenLineL = len(lineL)
        if header:
            if not col_index:
                try:
                    col_index = lineL.index(col_name)
                except ValueError:
                    print >>sys.stderr, "Un-exisited column name %s" % col_name
                    sys.exit(1)
            else:
                col_index -= 1
            #---------------------------------
            lineL.insert(col_index+1, "GeneNames")
            print >>file_out_fh, '\t'.join(lineL)
            header -= 1
            continue
        #-------------------------------------
        idL = lineL[col_index].split('/')
        nameL = '/'.join([idMapD.get(id) for id in idL])
        lineL.insert(col_index+1, nameL)
        print >>file_out_fh, '\t'.join(lineL)
    #-------------END reading file----------
    file_out_fh.close()
    os.system("/bin/mv -f "+file_out+" "+file)
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


