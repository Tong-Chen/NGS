#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2013, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Functional description:
    This is designed to transfer a normal expression data matrix to
    gct format file for GSEA usage.
    
    This will also generate a cls file.

    Ref: http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Expression_Data_Formats

    gct:
    1. The first line of gct is "#1.2"
    2. The second line of gct caontains two columns, the first column
    contains <the number of probe sets or genes in the file>, the
    second column contains <the number of samples >. Both should be
    less than the number of rows and columns in original file and new
    gct file.
    3. The third line is a header line with the first two columns
    containing <gene name> and <gene description (which is not used
    for GSEA)>. Third column onwards are unqiue sample names. 

    cls:
    1. The first line of cls is <(number of samples) (space) (number
    of classes) (space) 1>.
    2. The second line of cls is <# (space) (class 0 name) (space)
    (class 1 name) ... >.
    3. The third line of cls is <sample 1 class> <sample_2_class>


    Input file format:
    
    A.
    ==> test.A <==
    Gene	kdctcf_1	kdctcf_2	kdctcf_3	ctl_1	ctl_2	ctl_3
    a	1	2	3	4	5	6
    b	1	2	3	4	5	6
    c	1	2	3	4	5	6
    d	1	2	3	4	5	6

    RUN COMMAND: transferNormalExprMatrixForGSEA.py -i test.A

    ==> test.A.cls <==
    6 2 1
    # kdctcf ctl
    kdctcf kdctcf kdctcf ctl ctl ctl

    ==> test.A.gct <==
    #1.2
    4	6
    Gene	description	kdctcf_1	kdctcf_2	kdctcf_3	ctl_1	ctl_2	ctl_3
    a	a	1	2	3	4	5	6
    b	b	1	2	3	4	5	6
    c	c	1	2	3	4	5	6
    d	d	1	2	3	4	5	6


    B. (specify -u 0)
    ==> test.B <==
    Gene	kdctcf	kdctcf	kdctcf	ctl	ctl	ctl
    a	1	2	3	4	5	6
    b	1	2	3	4	5	6
    c	1	2	3	4	5	6
    d	1	2	3	4	5	6
    
    RUN COMMAND: transferNormalExprMatrixForGSEA.py -i test.B -u 0

    ==> test.B.cls <==
    6 2 1
    # kdctcf ctl
    kdctcf kdctcf kdctcf ctl ctl ctl

    ==> test.B.gct <==
    #1.2
    4	6
    Gene	description	kdctcf_1	kdctcf_2	kdctcf_3	ctl_4	ctl_5	ctl_6
    a	a	1	2	3	4	5	6
    b	b	1	2	3	4	5	6
    c	c	1	2	3	4	5	6
    d	d	1	2	3	4	5	6
    
    Both input files contain 6 samples and 2 classes with each class
    containing 3 samples. The string before the last underscore(_)
    will be treated as class names. 
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A normal expression matrix \
with the first line as header line contains sample names \
and first column contains gene names and other columns \
contain expression value. \
Pay attention to the header line. In default, strings before \
the last underscore(_) are used as class names.")
    parser.add_option("-u", "--unique-sample", dest="unique",
        metavar="unique", default=1, help="A variable to specify if \
each sample name is unique. Default 1 means each sample name is \
unique. If sample names have duplicate in given expression matrix \
file, please give <0> to <-u>, the program will add a unique label \
to each sample name to make them unique.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    prefix = file
    unique = int(options.unique)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    head = 1
    gct = open(prefix+'.gct', 'w')
    cls = open(prefix+'.cls', 'w')
    tmpL = []
    classL = []
    classLF = []
    row = 0
    col = 0
    for line in fh:
        if head:
            headL = line.strip().split('\t')
            len_headL = len(headL)
            col = len_headL - 1
            if not unique:
                for i in range(1, len_headL):
                    headL[i] = headL[i]+"_"+str(i)
            tmpL.append('\t'.join([headL[0], 'description','\t'.join(headL[1:])]))
            #get class information
            for i in range(1, len_headL):
                value = headL[i].rsplit('_', 1)[0]
                classLF.append(value)
                if value not in classL:
                    classL.append(value)
            head -= 1
            continue
        #----------------------------------
        lineL = line.strip().split('\t', 1)
        tmpL.append('\t'.join([lineL[0], lineL[0], lineL[1]]))
        row += 1
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    #-------BEGIN output-----------------
    print >>gct, "#1.2"
    print >>gct, "%d\t%d" % (row, col)
    print >>gct, '\n'.join(tmpL)
    gct.close()
    #------------------------------------
    print >>cls, "%d %d 1" % (col, len(classL))
    print >>cls, "# %s" % ' '.join(classL)
    print >>cls, ' '.join(classLF)
    cls.close()
    if verbose:
        print >>sys.stderr,\
            "--Successful %s" % strftime(timeformat, localtime())
if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()



