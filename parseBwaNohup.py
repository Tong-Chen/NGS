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
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #--------------------------------------------
    command = ('bwa aln', 'bwa samse')
    cmdDict = dict([(cmd,[]) for cmd in command])
    #runDict = {}
    runL = []
    errList = []
    killList = []
    for line in open('nohup.out'):
        if line[0] != '[':
            for cmd in command:
                if line.startswith(cmd):
                    file = (line.split('>')[-1]).strip()
                    cmdDict[cmd].append(file)
        elif line.startswith('[main]'):
                runL.append(line)
        elif line.startswith('make'):
            errList.append(line.strip())
        elif not line.startswith('[bwa_aln'):
            print line,
        if line.find('Killed') != -1:
            killList.append(line.strip())
        indexMake = line.find('make: ***')
        if indexMake != -1:
            errList.append(line[indexMake:-1])
    #------------------------------------------------
    print 'Startlist---------------------------------'
    for key, valueL in cmdDict.items():
        valueL.sort()
        print ">%s\n%s" % (key, '\t'.join(valueL))
    print 'Runlist---------------------------'
    print ''.join(runL)
    print 'Errlist----------------------------'
    print '\n'.join(errList)
    print 'killlist----------------------------'
    print '\n'.join(killList)
    #length = len(runL)
    #if length % 3:
    #    print "Uncomplete [main] line"
    #--------------------------------------
    #for i in range(1, length, 3):
        
if __name__ == '__main__':
    main()

