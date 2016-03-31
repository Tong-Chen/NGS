#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
#import sys
#sys.path.append("/home/ct/pylib")
'''


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
    import sys
    import re


    if len(sys.argv) != 2:
        print 'Wrong argument number.\nUsing python', sys.argv[0],'filename'
        sys.exit(0)
    else:
        large_symbol = re.compile('>(.*)')
        mail = re.compile('.+?@([0-9]+?)@.+?@([0-9]+?)@')
        adict = {}
        key = ' '
        poskey = ' '
        try:
            fh = open(sys.argv[1])
        except IOError, e:
            print sys.argv[0], 'open error', e
            sys.exit(0)
        else:
            #large_flag = 0
            #pos_flag = 0
            for line in fh:
                large = large_symbol.match(line)
                if large:
                    key = large.group(1)
                    adict[key] = {}
                pos = mail.search(line)
                if pos:
                    poskey, posvalue = pos.group(1), pos.group(2)
                    if poskey in adict[key]:
                        adict[key][poskey].append(posvalue)
                    else:
                        adict[key][poskey] = []
            fh.close()
            for key, value in adict.items():
                fh2 = open(key, 'w')
                poskeyArr = []
                for poskey, posvalue in value.items():
                    poskeyArr.append(poskey)
                    filename = key + '_' + poskey
                    fh = open(filename, 'w')
                    for eachposvalue in posvalue:
                        print >>fh, eachposvalue
                    fh.close()
                poskeyArr.sort()
                for eachposkey in poskeyArr:
                    print >>fh2, eachposkey
                fh2.close()


                



