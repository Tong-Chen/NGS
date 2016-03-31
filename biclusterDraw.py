#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
This is used to do ploting.

Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
==============================================================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = u'陈同(chentong) & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=============================================================================================
import matplotlib.pyplot as plt
import sys

def ct_plot(filename, fmt=5):
    #f = file(filename, 'r')
    #while True:
        #array =[int(i) for i in f.readline().split()]
	#array = [int(i) 
        #if len(array) == 0:
        #    break
        #plt.plot(array,'-.')
    #f.close()
    #plt.show()
    #fname = filename+'.png'
    #plt.savefig(fname)
    #plt.close()
    for line in open(filename):
        if line.find('-') != -1:
            fmt += 1
        array = [int(line[i:i+fmt]) for i in range(0, len(line)-1, fmt)]
        if len(array) == 0:
            break
        plt.plot(array,'-.')
    fname = filename+'.png'
    plt.savefig(fname)
    plt.close()


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "Need more argumnet. \nUsage: %s filename num of \
char in one column[positive integer, often the number from %5d]" % sys.argv[0]
        sys.exit()
    fh = file(sys.argv[1],'r')
    fmt = int(sys.argv[2])
    while True:
        filename = fh.readline().split()[0]
        if len(filename) < 0:
            break
        ct_plot(filename, fmt)
    fh.close()
