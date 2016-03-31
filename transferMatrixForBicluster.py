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
    print >>sys.stderr, "Transfer usual expression matrix to the \
format bicluster used, but with heads and labels in. This is the \
intermidate file which not only contains the numbers can be directly \
used in bicluster but also gene names to extract."
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 8:
        print >>sys.stderr, 'Using python %s input output head(0/num of \
head lines) label column(. seperated if more, 0 starts, \
only first several columns) output \
fmt(the number of characters in one column, similar with %s5f but only\
 number) time(an integer which related with the preceding parameter)\
 rep time[1 for no rep or you do not want to merge the rep.]' % (sys.argv[0], "%")
        sys.exit(0)
    #----------------------------------------
    output = sys.argv[2]
    fh = open(output, 'w')
    head = int(sys.argv[3])
    excluClo = [int(i) for i in sys.argv[4].split('.')]
    fmt = int(sys.argv[5])
    time = int(sys.argv[6])
    rep = int(sys.argv[7])
    first = 1
    for line in open(sys.argv[1]):
        if head:
            head -= 1
            print >>fh, line,
            #fh.write(line)
            continue
        #-------begin reading-------
        lineL = line.strip().split('\t')
        if first:  #compute the number of columns only the first time
            lenL = len(lineL)
        newLineL = []
        data = ''
        for i in excluClo:
            newLineL.append(lineL[i])
        for i in range(len(excluClo), lenL, rep):
            tmpLineL = [float(i52) for i52 in lineL[i:i+rep]]
            value = str(int(round(sum(tmpLineL) / rep * time, 0)))  
            lenValue = len(value)
            assert(fmt >= lenValue)
            data += ' ' * (fmt-lenValue) + value
        #----------------------------------------
        newLineL.append(data)
        print >>fh, '\t'.join(newLineL)
    #----------End reading----------------------------
    fh.close()
                
if __name__ == '__main__':
    main()

