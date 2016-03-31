#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
import sys
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
import urllib
def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s genelistfile others' % sys.argv[0]
#specie \
#condition format[if there are multiple species \
#(seperate with ":") or conditions, \
#or there are multiple words in one parameter, quote them ]' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    file = sys.argv[1]
    if len(file):
        geneList = [line.strip() for line in open(file)]
        if len(geneList):
            geneIs = 'geneIs=' + '+'.join(geneList)
        else:
            print >>sys.stderr, "Empty input gene list. If that is you\
want, using '' as file name."
            sys.exit()
    else:
        file = 'wholeGene'
        geneIs = ''
    #specie = sys.argv[2].split(':')
    others = sys.argv[2]
    if others[0] != '&':
        others = '&' + others
    #-----------------------------
    url =  "www.ebi.ac.uk/gxa/api/v1?" + geneIs + others
    url = url.replace('?&', '?', 1)
    print >>sys.stderr, url
    wp = urllib.urlopen(url)
    content = wp.read()
    fh = open(file+'.web', 'w')
    print >>fh, '#%s' % others
    fh.write(content)
    fh.close()
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


