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


def countCpG(file, substr):
    '''
    CG count
    >5 CpG per 100bp
    1-5 CpG per 100bp
    <1 CpG per 100bp
    '''
    fh = open('.'.join([file, '.'.join(substr), 'Num']),  'w')
    for line in open(file):
        if line[0] == '>':
            key = line[1:-1]
        else:
            line = line.strip()
            lenline = len(line)
            CpG = 0
            for i in substr:
                CpG += line.count(i)
            print >>fh, "%s\t%f" % (key, CpG * 1.0 / lenline)
        #===============================================
    #--------------------------------------------------------
    fh.close()
    #------------------------------------------------------
#------------------------------------------
def main():
    print >>sys.stderr, "Count the number of <substring> in the specified regions. \
If multiple substrings are given, it returns the number of all \
substrings related to one fasta sequence."
    print >>sys.stderr, "Print the result to file named after <file.substr1.substr2.Num>."
    if len(sys.argv) < 3:
        print >>sys.stderr, 'Using python %s fastafile substring1 [substring2] ' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------------------
    countCpG(sys.argv[1], sys.argv[2:])
if __name__ == '__main__':
    main()

