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
from ctIO import readRep
def main():
    print >>sys.stderr, "To detect the conservation among orthologs,\
use repetitions as the query and its related orthologs as db(after \
makeblastdb). "
    if len(sys.argv) != 3:
        print >>sys.stderr,'Using python %s repfile dbpath/' % sys.argv[0]
        sys.exit(0)
    #-------------------------------------
    file = sys.argv[1]
    if file.find('LCSs') != -1:
        label = '.LCSs'
    elif file.find('HCSs') != -1:
        label = '.HCSs'
    #patched at 20110922. Before not give the inital value to [label].
    #So it will give an error when dealing with non 'LCSs' and 'HCSs'
    #files.
    else:
        label = ''
    noOrtho = 0
    path = sys.argv[2]
    repDict = {}
    readRep(sys.argv[1], repDict)
    for locus, valueL in repDict.items():
        tmppath = path + locus
        #print tmppath
        if not os.path.exists(tmppath):
            noOrtho += 1
            continue
        #-------------------------------
        midlen = 30
        short = locus+label+'.short'
        long = locus +label+'.long'
        fhshort = open(short, 'w')
        fhlong = open(long, 'w')
        group = 0
        for groupD in valueL:
            group += 1
            tmpDict = {}
            groupDKeyL = groupD.keys()
            groupDKeyL.sort()
            for pos in groupDKeyL:
                seq = groupD[pos]
                if seq not in tmpDict:
                    tmpDict[seq] = [str(pos[0])]
                else:
                    tmpDict[seq].append(str(pos[0]))
            #-------------------------------------------
            tmpDictKeyL = tmpDict.keys()
            tmpDictKeyL.sort()
            for seq in tmpDictKeyL:
                lenseq = len(seq)
                pos = ':'.join(tmpDict[seq])
                if lenseq <= midlen:
                    print >>fhshort, '>%s.%s.%s\n%s' % \
                        (locus, str(group), pos, seq)
                else:
                    print >>fhlong, '>%s.%s.%s\n%s' % \
                        (locus, str(group), pos, seq)

            #--------END one group------------------------------
        fhshort.close()
        fhlong.close()
        cmdshort = ' '.join(('psiblast -query', short, '-db', tmppath, \
            '-out', short+'.out', '-num_iterations 5','-evalue 20000',\
            '-matrix PAM30', '-comp_based_stats 0', '-word_size 2'))
        cmdlong = ' '.join(('psiblast -query', long, '-db', tmppath, \
            '-out', long+'.out', '-num_iterations 5'))
        os.system(cmdshort)
        os.system(cmdlong)
        cmdshort = ' '.join(('psiblast -query', short, '-db', tmppath, \
            '-out', short+'.table', '-num_iterations 5','-evalue 20000',\
            '-matrix PAM30', '-comp_based_stats 0', '-word_size 2',
            '-outfmt 7'))
        cmdlong = ' '.join(('psiblast -query', long, '-db', tmppath, \
            '-out', long+'.table', '-num_iterations 5', '-outfmt 7'))
        #print cmd
        #break
        os.system(cmdshort)
        os.system(cmdlong)
        #------------END one locus
    print noOrtho
    #-------------END----all-----------------
if __name__ == '__main__':
    main()

