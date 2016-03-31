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
from optparse import OptionParser as OP

def cmdparameter(argv, options, args):
    usages = "%prog [-i blast] [-o output]"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="The result of blast default output")
    parser.add_option("-o", "--output-file", dest="fileout",
        metavar="FILEOUT", help="Save the parse")
    parser.add_option("-s", "--min-identify", dest="minI",
        metavar="MinIdentify", default="0.7", type="float",
        help="The min indentiy seted [%default]")
    parser.add_option("-l", "--max-identify", dest="maxI",
        metavar="MaxIdentify", default="1.0", type="float",
        help="The max indentiy seted [%default]")
    parser.add_option("-q", "--quiet", dest="verbose", default=True,
        help="Open the quiet mode[True]", action="store_false")
    parser.add_option("-d", "--debug", dest="verbose", default=False,
        help="Open the debug mode[%default]", action="store_true")
    #parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)
#--------------------------------------------------------------------

def main():
    options = {}
    args = []
    cmdparameter(sys.argv[1:], options, args)

if __name__ == '__main__':
    main()

