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
    print >>sys.stderr, "Output an eps and pdf graph file with the \
filename given as prefix."
    print >>sys.stderr, '''
File Format:(with a header line, startswith a '#'. This is necessary
and is used to annotate the first line not as data.)
#sample iso1    iso2    iso3    iso4
sample1 0.345   0.35    0.43    0.89
sample2 1.234   9.23    20.1    0
sample3 5.02    3.42    0.33    0.01
sample4 5.60    0.01    0.03    4.3
#######################
Note: sample column is the x-axis label
sample row is the legend
    '''
    if len(sys.argv) < 2 :
        print >>sys.stderr, 'Using python %s filename \
legend[default with legend, which is your itms in your first line, \
if there are many items, it would be better to tun it off by supply\
a 0 here.] xtics[default text xtics, if number supply a 0] \
preaditionallines[; seperate multiple lines] postaditionallines \
xlabel ylabel title [lines or linespoints]' % sys.argv[0]
        sys.exit(0)
    #---------------------------------
    lenpara = len(sys.argv)
    legendornot = 1
    if lenpara >= 3:
        legendornot = int(sys.argv[2])
    #-------------------------------
    xticsornot = 1 #text xtics
    xticsvalue = ''
    if lenpara >= 4:
        xticsornot = int(sys.argv[3])
        xticsvalue = '1:'
    #-------------------------------
    pregnu = ''
    if lenpara >= 5:
        pregnu = sys.argv[4]       
    #-------------------------------
    postgnu = ''
    if lenpara >= 6:
        postgnu = sys.argv[5]       
    #-------------------------------
    xlabel = ''
    if lenpara >= 7:
        xlabel = sys.argv[6]
    #-------------------------------
    ylabel = ''
    if lenpara >= 8:
        ylabel = sys.argv[7]
    #-------------------------------
    title = ''
    if lenpara >= 9:
        title = sys.argv[8]
    #-------------------------------
    style = 'lines'
    if lenpara >=10:
        style = sys.argv[9]
    #-------------------------------
    sample = []
    header = 1
    for line in open(sys.argv[1]):
        if header:
            header = 0
            lineL = line.strip().split()
            plot = 'plot ' + '\'' + sys.argv[1] + '\' '  
            pos = 1
            for item in lineL[1:]:
                pos += 1
                if pos == 2:
                    plot += 'using ' + xticsvalue + str(pos) + \
                        ' title ' + '\'' +\
                        item + '\'' + ' with ' + style
                else:
                    plot += ',\'\' using ' + xticsvalue + str(pos) + ' title ' + '\'' +\
                        item + '\'' + ' with ' + style
        #---------------------------------------
        else:
            tmp = line.split("\t",1)[0]
            sample.append(tmp)
    #--------------------------------------------
    xtics = 'set xtics ('
    i48 = -1
    for item in sample:
        i48 += 1
        if item == '-':
            continue
        if i48 == 0:
            xtics += '\'' + item + '\' ' + str(i48) 
        else:
            xtics += ',\'' + item + '\' ' + str(i48) 
    xtics += ')'
    #------------------------------
    xrange='set xrange [-0.5:' + str(i48+0.5) + ']'
    #------------------------------
    fileout = sys.argv[1] + '.plt'
    fh = open(fileout , 'w')
    print >>fh, 'set term postscript eps color'
    print >>fh, "set output \'" + sys.argv[1] + '.eps\''
    if not legendornot:
        print >>fh, "set nokey"
    if xticsornot:
        print >>fh, xrange
        print >>fh, xtics
    if pregnu:
        print >>fh, pregnu
    #----------------------------
    if xlabel:
        print >>fh, 'set xlabel \"%s\"' % xlabel 
    if ylabel:
        print >>fh, 'set ylabel \"%s\"' % ylabel 
    if title:
        print >>fh, 'set title \"%s\"' % title 
    #-----------------------------
    print >>fh, plot
    if postgnu:
        print >>fh, postgnu
    print >>fh, 'exit'
    fh.close()
    cmd1 = 'gnuplot ' + fileout
    cmd2 = 'convert -density 300 -flatten ' + sys.argv[1] + '.eps ' +\
        sys.argv[1] + '.png'
    os.system(cmd1)
    os.system(cmd2)
if __name__ == '__main__':
    main()

