#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
#import sys
#sys.path.append("/home/ct/pylib")
'''
this file is used to transfer "TAIR8_pep_20080412" for C use.

Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
==============================================================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=============================================================================================
if __name__ == '__main__':
    import re
    import sys
    print >>sys.stderr, 'this used to transfer at sequence for c'
    if len(sys.argv) < 2:
        print >>sys.stderr, 'Using %s filename' % sys.argv[0]
        sys.exit(1)

    #for there we use match(), which match at the beginning of the string. 
    begin_gt = re.compile('>(.{11,12}) |')
    #this means begin a new sequence
    flag = 0
    #fh2 = open("TAIR8_pep_20080412_FOR_C", 'w')
    for line in open(sys.argv[1]):
        if flag == 0:
            m = begin_gt.match(line)
            key = m.group(1)
            print key
            str = ''
            flag = 1
        else:
            if re.search('\*', line):
                flag = 0
                str += line[:-2]
                print len(str)
                print str
            else:
                str += line[:-1]
    #fh2.close()
