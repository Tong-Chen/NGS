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
    This is designed to get TCGA data for giving META manifest file from TCGA portal
    https://gdc-portal.nci.nih.gov/.

    1. Download data command like <gdc-client download -m TACG_GEQ_RPKM_UQ.gdc_manifest.2016-11-15T02_55_14.725068.tsv>
    2. Using GDC API to get sample information of UUIDs in META manifest file. See guide in <https://gdc-docs.nci.nih.gov/API/Users_Guide/>

META manifest file:

id	filename	md5	size	state
737a1837-5d14-45d7-af77-9b60fed2e556	fbd32fa7-a9c7-4148-a3d4-69179c833f97.FPKM-UQ.txt.gz	f4a1566f115f7ef6f6fd0e7cc2d9aebd	583252	submitted
fc10a611-a3e2-4ad6-8f76-970fa4f4f882	032c007c-a7c8-42bd-a7af-3134bcc8c4a1.FPKM-UQ.txt.gz	0382f24622fe77adce19130e726ad517	518849	submitted
7bb671f6-2569-481b-a136-e04789cce183	16c25e30-6e3d-401c-90c3-053af280fbb2.FPKM-UQ.txt.gz	a8607bc4aa598f1c799322251511929b	505732	submitted
eb44256d-5e56-4d8e-bbc8-d9e8055342d9	3490b791-16ec-41e8-9015-5f3d33c856c0.FPKM-UQ.txt.gz	882107dba4b4b277e7ee8b51f56e1aa5	510390	submitted
73f58881-2c4e-43cb-b69d-50e18d9cc909	7d6e7cb6-598a-4723-8ea1-60fee87edafe.FPKM-UQ.txt.gz	982f8ad48735fd7b9b1bf1ac0a6be39d	490914	submitted
629b971c-909b-4ccd-96f4-e53ef5a5ec7c	84aa891c-c73e-48eb-bdbe-5774c2dc3c58.FPKM-UQ.txt.gz	f9cb42c8f7f5308f53d2bf6e5276153b	513349	submitted
a966ac19-c902-46a1-a3a5-3be659ef569b	8ac283d9-5595-4378-8f4d-55c1d01ff9c0.FPKM-UQ.txt.gz	6b8527a1fca43ebfed7328f18f60d1e4	509802	submitted
ee8a0ff1-5267-49ee-ae74-ad2b2a04b36b	444e9b82-d8c5-41b8-a51b-45d9717fbce0.FPKM-UQ.txt.gz	4b831d83498f5743924522a8f625d1f7	512767	submitted
d7fddb99-6df0-4858-a563-7cbf3b99862a	c0b05393-e154-459d-a17a-9e03d404ae6e.FPKM-UQ.txt.gz	e574f0a1c7a8f0a35b51c5a69541660a	512913	submitted
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
        metavar="FILEIN", help="META manifest downloaded from GDC data portal")
    parser.add_option("-d", "--data-download", dest="data",
        default=False, action="store_true", help="Given <-d> to download data.")
    parser.add_option("-s", "--sample-download", dest="sample",
        default=False, action="store_true", help="Given <-s> to download sample info data.")
    parser.add_option("-v", "--verbose", dest="verbose",
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
    download_data = options.data
    download_sample = options.sample
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if download_data:
        cmd = "nohup gdc-client download -m " + file + ' &'
        os.system(cmd)
    if not download_sample:
        return

    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    header = 1
    idL = []
    for line in fh:
        if header:
            header -= 1
            continue
        id = line.strip().split('\t')[0]
        idL.append(id)
    #--------------------------------------
    idS = '"'+'","'.join(idL)+'"'
    #-------------END reading file----------
    payload = file+'.Payload.txt'
    payload_fh = open(payload, 'w')
    print >>payload_fh, '''{
    "filters": {
        "op": "in", 
        "content": {
            "field": "files.file_id", 
            "value": [%s]    
        }
    }, 
    "format": "TSV", 
    "fields": "file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id", 
    "size": "100"
}''' % idS
    payload_fh.close()
    cmd = 'curl --request POST --header "Content-Type: application/json" --data @{} "https://gdc-api.nci.nih.gov/files"  > {}.meta'.format(payload, file)
    if debug:
        print >>sys.stderr, cmd
    else:
        os.system(cmd)
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


