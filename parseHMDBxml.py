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
    This is designed to parse HMDB metabolites xml file.
    http://www.hmdb.ca/downloads/hmdb_metabolites.zip

    It will extract HMDB accession, cas_registry_number, 
    kegg_id, biocyc_id, 
    iupac_name, traditional_iupac,  synonyms,
    average_molecular_weight, chemical_formula, 
    smiles
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
from bs4 import BeautifulSoup
reload(sys)
sys.setdefaultencoding('utf8')
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
        metavar="FILEIN", help="HMDB_metabolites.xml. \
Only one compound are allowed in each xml file.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def null(string, null='NA'):
    if string == None:
        return null
    else:
        return string

#--------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    xml_soup = BeautifulSoup(open(file), 'xml')
    #for metabolism in xml_soup.find_all('metabolite'):
    hmdb_accession = null(xml_soup.accession.string)
    hmdb_name = null(xml_soup.find('name').string)
    cas_registry_number = null(xml_soup.cas_registry_number.string)
    chemspider_id = null(xml_soup.chemspider_id.string)
    kegg_id = null(xml_soup.kegg_id.string)
    biocyc_id = null(xml_soup.biocyc_id.string)
    iupac_name = null(xml_soup.iupac_name.string)
    traditional_iupac = null(xml_soup.traditional_iupac.string)
    average_molecular_weight = null(xml_soup.average_molecular_weight.string)
    chemical_formula = null(xml_soup.chemical_formula.string)
    smiles = null(xml_soup.smiles.string)
    synonyms = xml_soup.synonyms.strings
    synonymL = []
    if synonyms:
        for synonym in synonyms:
            #print synonym, 
            synonym = synonym.strip()
            #print synonym
            if synonym:
                #print 'here'
                synonymL.append(synonym)
        #-----------------------------------
    #------------------------------------
    #print synonymL
    print '\t'.join(['hmdb_accession', 'hmdb_name',
        'cas_registry_number', 'chemspider_id', 'kegg_id', 'biocyc_id', 'iupac_name',
        'traditional_iupac', 'chemical_formula', 'synonyms',
        'average_molecular_weight', 'smiles'])
    print '\t'.join([hmdb_accession, hmdb_name, 
        cas_registry_number, chemspider_id, kegg_id,
        biocyc_id, iupac_name, traditional_iupac,
        chemical_formula, '___'.join(synonymL),
        average_molecular_weight, smiles])
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


