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
    This is designed to classify TF families.

<family> file:
# Two pfam domains are separated by <;>
# <PF00847:2> represents at least two copies of PF00847 are required
#             for AP2 subfamily.
# <PF00847:1;PF02362:1> represents at least one copy of PF00847 and
#             PF02362 are required for RAV subfamily.
# <PF06943:1;PF00656:0> represents at least one copy of PF06943 and
#             without PF00656 are required for LSD subfamily.
Family	SubFamily	Standard
AP2/ERF	AP2	PF00847:2
AP2/ERF	RAV	PF00847:1;PF02362:1
BBR-BPC	BBR-BPC	PF06217:1
BES1	BES1	PF05687:1
C2C2	LSD	PF06943:1;PF00656:0
C3H	C3H	PF00642:1;PF00076:0;PF00271:0
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
        metavar="FILEIN", help="Pfam output or Trinotate annotation")
    parser.add_option("-t", "--input-file-type", dest="file_type",
        metavar="FILE-TYPE", help="pfam or trinotate")
    parser.add_option("-f", "--family", dest="family",
        help="TF family file as described above. \
For plant TFs, </MPATHB/resource/TFs/plant/PlantTFDB.family.std>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
#------------------------------------
def saveDict(aDict, key, value):
    if key not in aDict:
        aDict[key] = {}
    if value not in aDict[key]:
        aDict[key][value] = {}
    return aDict[key][value]
#------------------------------------

def readTFfamily(file):
    head = 1
    existDict = {}
    '''
    existDict = {
        'PF00643' : {1:{'PF06203':{1:{'C2C2    CO-like':'C2C2
        CO-like'}}}}
    }
    '''
    excludeDict = {}
    '''
    excludeDict = {
        'C2H2    C2H2': ['PF00929']
        'C3H     C3H' : ['PF00076', 'PF00271']
    }
    '''
    for line in open(file):
        if head:
            head -= 1
            continue
        #--------------------------
        lineL = line.strip().rsplit('\t', 1)
        family, pfam = lineL
        pfamL = pfam.split(';')
        pfamL.sort()

        pfam, count = pfamL[0].split(':')
        count = int(count)
        if count > 0:
            tmpD = saveDict(existDict, pfam, count)
        else:
            if family not in excludeDict:
                excludeDict[family] = []
            excludeDict[family].append(pfam)

        for pfam_count in pfamL[1:]:
            pfam, count = pfam_count.split(':')
            count = int(count)
            if count > 0:
                tmpD = saveDict(tmpD, pfam, count)
            else:
                if family not in excludeDict:
                    excludeDict[family] = []
                excludeDict[family].append(pfam)
        #------------------------------------
        tmpD[family] = family   
    #-----------------------------------------------
    return existDict, excludeDict
#-----END readTFfamily-----------------------

def readTrinotate(trinotate, name_col=0, pfam_col=11):
    geneD = {}
    '''
    geneD = {
        'gene1': {'PF000':[1, 'desp'], 'PF001':[2, 'desp']}, 
        'gene2': {'PF002':[1, 'desp'], 'PF003':[1, 'desp']}, 
    }
    '''
    head = 1
    if head:
        head -= 1
        continue
    #--------------------------
    lineL  = line.split('\t')
    gene   = lineL[name_col]
    assert gene not in geneD, "Duplicate %s" % gene
    geneD[gene] = {}
    pfam_c = lineL[pfam_col].split('`')
    for pfam_i in pfam_c:
        pfam，desp1, desp2 = pfam_i.split('^')
        pfam = pfam.split('.')[0]
        if pfam not in geneD[gene]:
            geneD[gene][pfam] = [1, desp1]
        else:
            geneD[gene][pfam][0] = geneD[gene][pfam][0] + 1
    return geneD
#-------------------------------------

def pfamFamily(geneD, existDict, excludeDict):
    familyD = {}
    for gene, pfamD in geneD.items():
        pfamK = pfamD.keys()
        pfamK.sort()
        for pfam in pfamK:
            if pfam not in existDict:
                
    #---------------------------------------
#--------------pfamFamily----------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    anno = options.filein
    type = options.file_type
    if type == 'trinotate':
        geneD = readTrinotate(anno)
    elif type == 'pfam':
        geneD = readPfam(anno)
    else:
        print >>sys.stderr, "Unknown annotation type"
    #-----------------------------------------
    family = options.family
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    existDict, excludeDict = readTFfamily(family)
    #print existDict
    #print excludeDict
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


