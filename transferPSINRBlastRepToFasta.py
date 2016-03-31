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
    print >>sys.stderr, "Print the parts result to screen"
    print >>sys.stderr, "Transfer the psiblast hitted results to fasta\
and construct hmm."
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s filename outpath/' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------------------------
    fh = ''
    path = sys.argv[2]
    for line in open(sys.argv[1]):
        if line[0] == '=':
            if fh:
                fh.close()
            file = line[1:].split()[0]
            #print file
            tmpfile = file
            file = path + file
            fh = open(file, 'w')
        elif line[0] == '>':
            locus = line[1:-1]
        else:
            lineL = line.strip().split('#')
            for item in lineL:
                seq, pos = item.split(':')
                print >>fh, ">%s_%s\n%s" % (locus, pos, seq)
                print ">%s_%s_%s\n%s" % (tmpfile, locus, pos, seq)
    #---------------END reading------------------
    if fh:
        fh.close()
    #--------------The last item---------------------

if __name__ == '__main__':
    main()

