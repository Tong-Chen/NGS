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
import os

def main():
    print >>sys.stderr, "Download data from ncbi or ebi with ascp"
    lensa = len(sys.argv)
    if lensa < 2:
        print >>sys.stderr, 'Using python %s url[point to a file or \
directory] type[default ebi. \
another choice is ncbi, ori] \
dest[default ./]' % sys.argv[0]
        sys.exit(0)
    #------------------------------------
    url=sys.argv[1]
    if lensa == 3:
        type = sys.argv[2]
    else:
        type = 'ebi'
    if len(sys.argv) == 4:
        dest = sys.argv[3]
    else:
        dest = '.'
    #-------------------------------
    if not url.startswith("ftp://ftp"):
        print >>sys.stderr, "Wrong format, should be ftp://ftp.sra.ebi.ac.uk/"
        print >>sys.stderr, "Wrong URL:", url
    if type == 'ebi':
        url = url.replace('ftp://ftp', 'era-fasp@fasp')
        url = url.replace('/', ':/', 1)
    elif type == 'ncbi':
        url = url.replace('ftp://ftp-trace.ncbi.nlm.nih.gov', \
            'anonftp@ftp-private.ncbi.nlm.nih.gov:')
        if url.find('anonftp@') == -1:
            print >>sys.stderr, "Wrong urls. Transfer by your self. \
                Substitute things before the first single / \
                with anonftp@ftp-private.ncbi.nlm.nih.gov:"
            sys.exit(1)
    cmd = \
        'ascp -QTr -k1 -l 200M -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty '\
        + url + ' ' + dest + ' | tee >ascp.log 2>&1'
    print cmd
    os.system(cmd)
if __name__ == '__main__':
    main()

