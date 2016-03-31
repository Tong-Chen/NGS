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

def main():
    if len(sys.argv) < 3:
        print >>sys.stderr, "Extract sequences from fastq files.\
It was first planned for extract unmapped sequences \
and delete adaptor. Besides, one can also specify to \
output fasta files."
        print >>sys.stderr, "Print the result to screen"
        print >>sys.stderr, 'Using python %s fastq name out_fasta' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------------------------------
    aDict = {}
    for name in open(sys.argv[2]):
        aDict[name.strip()] = 1
    #--------------------------------
    if len(sys.argv) == 4:
        out_fasta = 1
    else:
        out_fasta = 0
    #----------------------
    output = 0
    i = 1
    for line in open(sys.argv[1]):
        if i % 4 == 1:
            if len(aDict) == 0:
                break
            key=line.split()[0].replace('@', '', 1)
            if key in aDict:
                output = 1
                aDict.pop(key)
            #else:
            #    print >>sys.stderr, key
        elif i % 4 == 0:
            if output and not out_fasta:
                print line,
            output = 0
        if output:
            if out_fasta:
                if i % 4 == 1:
                    print ">%s" % line,
                elif i % 4 == 2:
                    print line,
            else:
                print line, 
        i += 1
    #----------------------------------------------------------
#-------------------------------------------------------------------
if __name__ == '__main__':
    main()

