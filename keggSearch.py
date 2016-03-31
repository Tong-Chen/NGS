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

    This is designed to use togows to access KEGG database. 

    Given a <pathway> ID, one can get <orthologs> in this pathway,
    <genes> in each <orthology>, <aafasta> (amino acid sequences) of
    these <genes>.

    Given a <pathway> ID, one can get <compounds> in this pathway, 
    <names> (synonym) of each <compound>.

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime, sleep 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from bs4 import BeautifulSoup
import requests

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
        metavar="FILEIN", help="A file containing a list of IDs. \
like <ko00010> or <K00001> or <dme:Dmel_CG3481> or <C00003>")
    parser.add_option("-t", "--input-type", dest="in_type",
        help="The type of IDs in file given to <-i>, currently \
[pathway, orthology, genes, compound ] are supported.")
    parser.add_option("-T", "--out-type", dest="out_type",
        help="The type outputs, currently \
[orthologs, genes, aafasta, compounds, names, compound_formula_mw] are supported.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def getJson(in_type, out_type, id, try_time=5):
    url = ["http://togows.org/entry/kegg-", in_type, "/", id, "/",
        out_type, ".json"]
    #print >>sys.stderr, ''.join(url)
    data = requests.get(''.join(url))
    if data.status_code == requests.codes.ok:
        return data.json()
    else:
        print >>sys.stderr, ''.join(url)
        return ['NA']

#-----------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file     = options.filein
    in_type  = options.in_type
    out_type = options.out_type
    verbose  = options.verbose
    global debug
    debug    = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    for line in fh:
        id = line.strip()
        print >>sys.stderr, "Working on %s" % id
        if in_type == "pathway":
            if out_type in ["orthologs", "compounds"]:
                orthologD = getJson(in_type, out_type, id)[0]
                for k, d in orthologD.items():
                    print "%s\t%s" % (k, d)
            elif out_type == "genes":
                orthologD = getJson(in_type, "orthologs", id)[0]
                for k, d in orthologD.items():
                    geneD = getJson("orthology", out_type, k)[0]
                    for species, nameL in geneD.items():
                        for name in nameL:
                            print "\t".join([id, k, d, species, name])
            elif out_type == "aafasta":
                orthologD = getJson(in_type, "orthologs", id)[0]
                for k, d in orthologD.items():
                    geneD = getJson("orthology", "genes", k)[0]
                    geneIDL = []
                    for species, nameL in geneD.items():
                        for name in nameL:
                            print "\t".join([id, k, d, species, name])
                        geneIDL.extend([':'.join([species, name]) for name in nameL])
                    #geneIDs = ','.join(geneIDL)
                    output = open(id + '.' + k + '.fa', "w")
                    len_geneIDL = len(geneIDL)
                    for i in range(0, len_geneIDL, 10):
                        geneIDs = ','.join(geneIDL[i:i+10])
                        fastaL = getJson("genes", out_type, geneIDs)
                        print >>output, '\n'.join(fastaL)
                    output.close()
                    sleep(60)
            elif out_type == "names":
                orthologD = getJson(in_type, "compounds", id)[0]
                for k, d in orthologD.items():
                    compoundL = getJson("compound", out_type, k)[0]
                    for compound in compoundL:
                        print "\t".join(id, k, d, compound)
            elif out_type == "compound_formula_mw":
                outputL = getJson(in_type, "compounds", id)
                if len(outputL) == 0:
                    print >>sys.stderr, "No requests for %s" % ' '.join([in_type,
                        'compounds', id])
                    sleep(10)
                    continue
                orthologD = outputL[0]
                count = 0
                if orthologD:
                    for k, d in orthologD.items():
                        formula = getJson("compound", "formula", k)[0]
                        exact_mass = str(getJson("compound", "exact_mass", k)[0])
                        print "%s\t%s\t%s\t%s" % (k, d, formula, exact_mass)
                        count += 1
                        if count % 20 == 0:
                            sleep(count)
                    #------------END tracing orthologD--------
                else:
                    print >>sys.stderr, "No requests for %s" % id
                #---------Judge orthologD---------------------
                if count < 60: count = 60
                if count > 180: count = 180
                sleep(count)
            else:
                print >>sys.stderr, "Unsupported %s" % out_type
            #-------------------------------------
        elif in_type == "orthology":
            if out_type == "genes":
                geneD = getJson(in_type, out_type, id)[0]
                for species, nameL in geneD.items():
                    for name in nameL:
                        print "\t".join([id, species, name])
            elif out_type == "aafasta":
                geneD = getJson(in_type, "genes", id)[0]
                geneIDL = []
                for species, nameL in geneD.items():
                    for name in nameL:
                        print "\t".join([id, species, name])
                    geneIDL.extend([':'.join([species, name]) for name in nameL])
                output = open(id + '.fa', "w")
                len_geneIDL = len(geneIDL)
                for i in range(0, len_geneIDL, 10):
                    geneIDs = ','.join(geneIDL[i:i+10])
                    fastaL = getJson("genes", out_type, geneIDs)
                    print >>output, '\n'.join(fastaL)
                output.close()
            else:
                print >>sys.stderr, "Unsupported %s" % out_type
            #---------------------------------------------
        elif in_type == "genes":
            if out_type == "aafasta":
                fastaL = getJson("genes", out_type, id)
                print '\n'.join(fastaL)
            else:
                print >>sys.stderr, "Unsupported %s" % out_type
        elif in_type == "compound":
            if out_type == "names":
                orthologD = getJson(in_type, "compounds", id)[0]
                for k, d in orthologD.items():
                    compoundL = getJson("compound", out_type, k)[0]
                    for compound in compoundL:
                        print "\t".join(id, k, d, compound)
            elif out_type == "compound_formula_mw":
                formula = getJson("compound", "formula", id)[0]
                exact_mass = str(getJson("compound", "exact_mass", id)[0])
                print "%s\t%s\t%s" % (id, formula, exact_mass)
            else:
                print >>sys.stderr, "Unsupported %s" % out_type
            #-------------------------------------
        else:
            print >>sys.stderr, "Unsupported %s" % in_type
        #-------------------------------------------------------

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


