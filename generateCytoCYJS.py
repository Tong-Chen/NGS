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
    This is designed to generate cytoscape JSON file by supplying a three column network file and a node attribute file.

network file:

samp1	samp2	pcor
A	B	0.5
A	C	0.6
A	D	0.7
A	E	-0.2
A	F	-0.3
A	G	-0.6
A	H	-0.8
A	I	1

node attribute file

1. Currently,  only <SELF>, <POS>, <NEG> are supported. 
   <SELF> will be put at the center line.
   <POS> will be put at the top half circle.
   <NEG> will be put at the bottom half circle.

ID	Type
A	SELF
B	POS
C	POS
D	POS
E	NEG
F	NEG
G	NEG
H	NEG
I	POS

'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
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
        metavar="FILEIN", help="Network file with format specified above.")
    parser.add_option("-n", "--node", dest="node_attr",
        help="Node attribute file with format specified above.")
    parser.add_option("-o", "--output", dest="output",
        help="Output filename")
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
    attr = options.node_attr
    output = options.output
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------

    aDict = {}
    shared_name = file
    name = file

    aDict['format_version'] = "1.0"
    aDict["generated_by"] = "cytoscape-3.5.0"
    aDict["target_cytoscapejs_version"] = "~2.1"
    aDict["data"] = {}
    aDict["data"]["shared_name"] = "a.txt"
    aDict["data"]["name"] = "a.txt"
    aDict["data"]["SUID"] = 62
    aDict["data"]["__Annotations"] = []
    aDict["data"]["selected"] = 'true'
    aDict["elements"] = {}
    aDict["elements"]["nodes"] = []

    header = 1
    id = 0
    idD = {}
    for line in open(attr):
        lineL = line.strip().split('\t')
        if header:
            headerL = lineL[:]
            header -= 1
            continue
        subD = {}
        subD["data"] = {}
        id += 1
        name = lineL[0]
        idD[name] = id
        subD["data"]["id"] = str(id)
        subD["data"]["SUID"] = id
        subD["data"]["shared_name"] = name
        subD["data"]["name"] = name
        subD["data"]["selected"] = "false"
        for type, value in zip(headerL[1:], lineL[1:]):
            subD["data"][type] = value
        subD["position"] = {"x":0, "y":0}
        subD["selected"] = "false"
        aDict["elements"]["nodes"].append(subD)
    #----------------------------------------
    aDict["elements"]["edges"] = []
    
    header = 1
    for line in fh:
        lineL = line.strip().split('\t')
        if header:
            headerL = lineL[:]
            header -= 1
            continue
        subD = {}
        subD["data"] = {}
        source, target, pcor = lineL
        pcor = float(pcor)
        id += 1
        subD["data"]["id"] = str(id)
        subD["data"]["SUID"] = id
        subD["data"]["source"] = str(idD[source])
        subD["data"]["target"] = str(idD[target])
        name = " (interacts with) ".join(lineL[:2])
        subD["data"]["shared_name"] = name
        subD["data"]["name"] = name
        interaction = "interacts with"
        subD["data"]["shared_interaction"] = interaction
        subD["data"]["interaction"] = interaction
        subD["data"]["selected"] = "false"
        subD["data"][headerL[2]] = pcor
        subD["selected"] = "false"
        aDict["elements"]["edges"].append(subD)


    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    output_fh = open(output, 'w')
    print json_dump(aDict, output_fh, indent=4)
    output_fh.close()
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


