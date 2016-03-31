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
import pylab
from Bio import Phylo

def main():
    print >>sys.stderr, "Print the result to file"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s file.tre[nexus]' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------------------
    file = sys.argv[1]
    progN = 'twopi'
    #progN = 'neato'
    tree = Phylo.read(file, 'nexus')
    Phylo.draw_graphviz(tree, prog=progN)
    file2 = file.replace('tre', 'png')
    pylab.savefig(file2)
if __name__ == '__main__':
    main()

