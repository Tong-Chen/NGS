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
    This is designed to fetch gb file using entrez to get CDS and
    protein sequences. 
    This will replace entrezFetch.py
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime, sleep 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from Bio import Entrez
Entrez.email = "chentong_biology@163.com"

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
        metavar="FILEIN", help="A file containing a \
list of accessions (one in each line) or GI numbers")
    parser.add_option("-d", "--database", dest="db",
        default="nuccore", help="Default `nuccore`.")
    parser.add_option("-o", "--output-prefix", dest="out_pre",
        help="Prefix for output files (both protein and \
CDS sequences will be output.)")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    db = options.db
    out_pre = options.out_pre
    verbose = options.verbose
    global debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    
    prot_fh = open(out_pre+".prot.fa", 'w')
    nucl_fh = open(out_pre+".nucl.fa", 'w')

    count = 0
    for line in fh:
        acc = line.strip()
        count += 1
        sleep(count / 20)
        gb = Entrez.efetch(db=db, rettype="gb", id=acc, retmode='xml')
        records = Entrez.parse(gb)
        records = list(records)
        record  = records[0]
        dna     = record['GBSeq_sequence']
        prot_id = 'None'
        accession_ver = record['GBSeq_accession-version']
        GBSeq_feature_table = records[0]['GBSeq_feature-table']
        cds = prot = ''
        for feature in GBSeq_feature_table:
            if feature['GBFeature_key'] == 'mRNA':
                GBFeature_quals = feature['GBFeature_quals']
                transcription = 0
                for GB_qualsD in GBFeature_quals:
                    if GB_qualsD['GBQualifier_name'] == 'transcription':
                        transcription = 1
                        cds = GB_qualsD['GBQualifier_value']
                        break
                assert transcription, "No CDS found for %s" % acc
            elif feature['GBFeature_key'] == 'CDS':
                GBFeature_quals = feature['GBFeature_quals']
                GBFeature_location = feature['GBFeature_location']
                GBFeature_location = GBFeature_location.replace('join(', '')
                GBFeature_location = GBFeature_location.replace(')', '')
                GBFeature_locationL = GBFeature_location.split(',')
                cdsL = []
                for GBFeature_location in GBFeature_locationL:
                    start, end = GBFeature_location.split('..')
                    start = start.replace('<', '')
                    end   = end.replace('>', '')
                    cdsL.append(dna[(int(start)-1):int(end)])
                cds2 = ''.join(cdsL).upper()
                if cds:
                    try:
                        assert cds2 in cds, acc
                    except AssertionError:
                        print >>sys.stderr, "CDS direct:", cds
                        print >>sys.stderr, "CDS extract:", cds2
                        print >>sys.stderr, acc
                #---------------------------------------------
                cds = cds2
                translation = 0
                for GB_qualsD in GBFeature_quals:
                    if GB_qualsD['GBQualifier_name'] == 'translation':
                        translation = 1
                        prot = GB_qualsD['GBQualifier_value']
                    elif GB_qualsD['GBQualifier_name'] == 'db_xref':
                        gi = GB_qualsD['GBQualifier_value']
                    elif GB_qualsD['GBQualifier_name'] == 'protein_id':
                        prot_id = GB_qualsD['GBQualifier_value']
                
                assert translation, "No protein found for %s" % acc
            #------------------------------
        #---------------------------------------------
        print >>prot_fh, ">%s\t%s\t%s\n%s" % (acc, accession_ver, prot_id, prot)
        print >>nucl_fh, ">%s\t%s\t%s\n%s" % (acc, accession_ver, prot_id, cds)

    #-------------END reading file----------
    prot_fh.close()
    nucl_fh.close()
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


