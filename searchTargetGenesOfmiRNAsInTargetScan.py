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
    This is designed to extract target genes of specific miRNAs
    predicted by TargetScan. The full table and alignment file will be
    outputted.
'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import urllib2
import re
debug = 0

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
        metavar="FILEIN", help="The file contains miRNA names \
with each at one row.")
    parser.add_option("-s", "--species", dest="spec",
        default='Human', help="Specify the species information, \
default Human, accept Rat, Mouse or other supported species.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def parseContent(contentL, mir):
    pattern = re.compile('(<[^>]*>)*([^><]*)')
    next_pat = re.compile('"(.*)"')
    tmpL = []
    tmpL.append([])
    tmpL.append([])
    begin = 0
    for line in contentL:
        if line.find('Target gene') != -1:
            begin = 1
        #------------------------------
        if begin and line.startswith('<td>'):
            if debug:
                print >>sys.stderr, "### Get to expected miRNA-gene table"
            matchAll = pattern.findall(line)
            gene = matchAll[0][1]
            tr   = matchAll[1][1]
            func = matchAll[2][1]
            #mir  = matchAll[-6][1]
            tmpL[0].append(gene+'\t'+mir)
            if debug:
                print >>sys.stderr, '###', line
                print >>sys.stderr, '###', gene, mir
            next_url = next_pat.search(matchAll[-3][0]).groups()[0]
            alignL = getAlignment(next_url, mir)
            tmpL[1].append(alignL)
        #-------END one target gene-----------
    #-----------END all genes-----------------
    if begin == 0:
        print >>sys.stderr, "No miRNA %s in databse" % mir
    return tmpL
#---------------------------------------------------------

def getAlignment(next_url, mir):
    alignL = []
    next_contentL = urllib2.urlopen(next_url).read().split('\n')
    if debug:
        print >>sys.stderr, '\n'.join(next_contentL)
    pattern = re.compile('<.*?>([^><]+)(<.*?>)*([^><]+)')
    pattern_no = re.compile('<[^>]*>') #the substituted strings
    begin = 0
    for i in next_contentL:
        if i.startswith('<TD nowrap>') and i.find(mir) != -1:
            begin = 1
            if debug:
                print >>sys.stderr, '###Aln', i, 
            gene, deplete, mir_2 = pattern.search(i).groups()
            len_gene = len(gene)
            len_mir_2 = len(mir_2)
            expect_len = len_gene if len_gene > len_mir_2 else len_mir_2
            if debug:
                print >>sys.stderr, '###Aln', gene, mir_2
            
            if mir_2 != mir:
                print "### Unconsistent miRNA: \
original %s, here %s" % (mir, mir_2)
            #assert mir_2 == mir, "Unconsistent miRNA: \
#original %s, here %s" % (mir, mir_2)
        elif begin:
            assert i.startswith('<TD nowrap>')
            begin = 0
            pre, mid, post, unuse = i.split('<br>')
            pre = pattern_no.sub('', pre).replace('&nbsp;', ' ')
            mid = pattern_no.sub('', mid).replace('&nbsp;', ' ')
            post = pattern_no.sub('', post).replace('&nbsp;', ' ')
            if debug:
                print >>sys.stderr, '###', i
                print >>sys.stderr, '###', pre, mid, post
            alignL.append([gene+'\t'+pre, ' ' * expect_len + '\t'+mid, mir_2+ ' '*
                    (expect_len-len_mir_2) + '\t'+post])
        #-----------------------------------------------
    return alignL
#----------------------------------------------------------

def output(stringL):
    '''
    stringL = [ [gene\tmir, 
                 gene2\tmir
                ], 
                [ [
                    [gene\tseq, mir\tseq], 
                    [gene\tseq, mir\tseq], 
                  ], 
                  [
                    [gene2\tseq, mir\tseq],
                    [gene2\tseq, mir\tseq],
                  ]
                ]
              ]
    '''
    print '\n'.join(stringL[0])
    for listL in stringL[1]:
        for list in listL:
            print '\n'.join(list)
        print


#------------------------------------------------------------
def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    spec = options.spec
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    url_pre='http://www.targetscan.org/cgi-bin/targetscan/vert_61/targetscan.cgi?species='+ spec + '&gid=&mir_sc=&mir_c=&mir_nc=&mirg='
    for line in fh:
        mirna = line.strip()
        url = url_pre + mirna
        content = urllib2.urlopen(url).read()
        if debug:
            print >>sys.stderr, content 
        output(parseContent(content.split('\n'), mirna))
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
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



