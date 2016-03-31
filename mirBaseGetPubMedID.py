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
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
import urllib
import re
def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s url' % sys.argv[0]
        sys.exit(0)
    #-----------------------------------
    url = sys.argv[1]
    #---------get start-------------------
    urlL = url.split('&')
    lenurlL = len(urlL)
    for i in range(lenurlL):
        item = urlL[i]
        if item.find('start') == 0:
            start = int(item.split('=')[1])
            #startPos = i
            break
    #------get start---------------------------
    print >>sys.stderr, url
    wp = urllib.urlopen(url)
    content = wp.read()
    pat = re.compile('\["MGI:.+?"\]')
    mgi = pat.findall(content)
    currentGet = len(mgi)
    for item in mgi:
        print item[2:-2]
    while currentGet:
        oldstart = 'start='+str(start)
        start += currentGet -1 
        assert(start != 0)
        newstart = 'start='+str(start)
        assert(newstart != oldstart)
        url = url.replace(oldstart, newstart, 1)
        print >>sys.stderr, url
        try:
            wp = urllib.urlopen(url)
        except IOError:
            wp = urllib.urlopen(url)
        content = wp.read()
        mgi = pat.findall(content)
        currentGet = len(mgi)
        for item in mgi:
            print item[2:-2]
        #-------------------------------
    #--------------------------------------------
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()


