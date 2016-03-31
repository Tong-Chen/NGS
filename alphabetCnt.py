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


def countCpG(file, countAlpha):
    '''
    '''
    fh = open(file + '.' + '-'.join(countAlpha),  'w')
    print >>fh, "%s\t%s\t%s" % (file, '\t'.join(countAlpha),"total_len")
    for line in open(file):
        if line[0] == '>':
            key = line[1:-1]
        elif line[0] != '>':
            lenline = len(line) - 1
            cntL = [str(line.count(item)) for item in countAlpha]
            print >>fh, "%s\t%s\t%d" % (key, '\t'.join(cntL), lenline)
        #===============================================
    #--------------------------------------------------------
    fh.close()
#------------------------------------------
def main():
    if len(sys.argv) < 3:
        print >>sys.stderr, "Count the number of given alphabet(s) in \
the specified regions. Multiple alphabet(s) should be separated by -"
        print >>sys.stderr, 'Using python %s seq1-seq2-seq3.. \
fasta1[, fasta2...]' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------------------
    countAlpha = sys.argv[1].split('-')
    bedL = sys.argv[2:]
    for bed in bedL:
        countCpG(bed, countAlpha)
if __name__ == '__main__':
    main()

