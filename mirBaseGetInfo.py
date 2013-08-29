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
Functionla description

This is first designed to get the PUBMED ID from mirBase webpage.

'''
import urllib
import re
import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

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
        metavar="FILEIN", help="The file contains miRNA mature names \
(mirBase offical symbol).")
    parser.add_option("-o", "--operation", dest="op",
        default="pubmedID", metavar="pubmedID", help="Accept a \
string to represent the things you want. Default pubmedID. \
Currently, only this is accepted.")
    parser.add_option("-s", "--sleep-seconds", dest="sleep_s",
        default=1, help="The interval time between fetching. \
Sometimes this should be larger to avoid to be treated as attack. \
Default 1 second.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    op = options.op
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    url_p = "http://www.mirbase.org/cgi-bin/mature.pl?mature_acc="
    pat = re.compile(r"PMID:<a href=\"http://www.ncbi.nlm.nih.gov/pubmed/([0-9]+)")
    for line in fh:
        name = line.strip()
        url = url_p + name
        wp = urllib.urlopen(url)
        content = wp.read()
        valueL = pat.findall(content)
        print name, "\t".join(valueL)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
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

