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
    print >>sys.stderr, "------Print the result to file.------"
    if len(sys.argv) != 4:
        print >>sys.stderr, 'Using python %s groupfile speciesname \
acrony' % sys.argv[0]
        print >>sys.stderr, 'python %s groupfile \
Dicot.Arabidopsis.thaliana D.A.th' % sys.argv[0]
        sys.exit(0)
    
    count = 0
    dictLine = {}
    dictSpecies = {}
    speciesname = sys.argv[2]

    for line in open(sys.argv[1]):
        if line[0] == '>':
            count += 1
            dictLine[count] = ''
            dictSpecies[count] = set()
            continue
        elif line[0].isdigit() and line.find('|') != -1:
            species = line.split()[3]
            if species == speciesname:
                dictSpecies[count].add(speciesname)
        dictLine[count] += line
    #------------------------------------------------
    count2 = 0
    output = sys.argv[1] + '.' + sys.argv[3]
    fh = open(output, 'w')
    oldhandle = sys.stdout
    sys.stdout = fh
    for key, value in dictSpecies.items():
        if len(value) > 0:
            if 1:
                print '>' 
                print dictLine[key],
            count2 += 1
    #-----------------------------------------------
    sys.stdout = oldhandle
    fh.close()
    print >>sys.stderr, "There are %d groups original, now thre are\
 %d groups contain %s ." % (count, count2, speciesname)

if __name__ == '__main__':
    main()

