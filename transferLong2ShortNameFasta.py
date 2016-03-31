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

def labelP(num, alphabet):
    division = num / 26
    remain = num % 26
    label = ''
    label = alphabet[remain] + label
    while division > 25:
        remain = division % 26
        division = division / 26
        label = alphabet[remain-1] + label
    #---------------------------
    if division > 0:
        label = alphabet[division-1] + label
    return label
#-----------------------------------


def main():
    print >>sys.stderr, "Print the result to screen, and a map file \
with_a suffix .Name_Map."
    print >>sys.stderr, "Transfer sequence locus from long name to\
short name for fasta format file."
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------------
    alphabet = [chr(i) for i in range(97,123)]
    i = 0
    map40 = sys.argv[1] + '.Name_Map' 
    fh = open(map40, 'w')
    for line in open(sys.argv[1]):
        if line[0] == '>':
            locus = line[1:-1]
            newlocus = labelP(i, alphabet)
            print '>%s' % newlocus
            print >>fh, '%s\t%s' % (locus, newlocus)
            i += 1
        else:
            print line,
    fh.close()
if __name__ == '__main__':
    main()

