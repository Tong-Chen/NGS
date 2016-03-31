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
from ctIO import readRep

def readInterpro(interpro):
    '''
    aDict = {
        locus : [(pos1,pos2), (pos1,pos2)]
    }
    '''
    aDict = {}
    for line in open(interpro):
        lineL = line.split('\t',8)
        locus = lineL[0]
        if locus not in aDict:
            aDict[locus] = []
        aDict[locus].append((int(lineL[6]), int(lineL[7])))
    #-----------------------------------------
    return aDict

def output(repDict):
    '''
    repDict = 
    {
    AT1G62760: 
        [
            {(25, 31): 'SSLSPSS:25!', (51,57): 'SSLSPSSi:51!'},
            {(52, 60): 'SLSPSSPPPi:52@',},    
        ] 
    }
    '''
    locusL = repDict.keys()
    locusL.sort()
    for locus in locusL:
        valueL = repDict[locus]
        print '>%s' % locus
        for dictrep in valueL:
            repv = dictrep.values()
            repv.sort(key=lambda s: int(s.split(':')[1][:-1]))
            print '#'.join(repv)
#-------------------------------------------------

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) != 3:
        print >>sys.stderr, 'Using python %s repname interpro'\
            % sys.argv[0]
        sys.exit(0)
    #---------------------------------------------
    aDict = readInterpro(sys.argv[2])
    repDict = {}
    '''
    repDict = 
    {
    AT1G62760: 
        [
            {(25, 31): 'SSLSPSS', (51,57): 'SSLSPSS'},
            {(52, 60): 'SLSPSSPPP',},    
        ] 
    }
    '''
    readRep(sys.argv[1], repDict)
    if 0:
        print aDict
        print repDict
    '''
    Some symbol:
    !: new
    ------      domain
           ---  rep
    @: overlap
    ------    domain
      ------- rep
    ^: rep in domain
    -----------    domain
      -------      rep
    $: domain in rep
       -------     domain
    -------------  rep
    '''


    for locus, valueL in repDict.items():
        if locus in aDict:
            hasDomain = 1
            domainPosL = aDict[locus]
        else:
            hasDomain = 0
        for dictrep in valueL:
            for posset in dictrep.keys():
                begin = posset[0]
                end = posset[1]
                if 0:
                    print "begin is %d, end is %d" % (begin, end)
                if hasDomain:
                    nooverlap = 1
                    for domainset in domainPosL:
                        ds = domainset[0]
                        de = domainset[1]
                        if 0:
                            print "ds is %d, de is %d" % (ds, de)
                        if (begin > ds and begin < de and end > de)\
                           or (begin < ds and end > ds and end < de):
                            dictrep[posset] += ':'+str(begin)+'@'
                            nooverlap = 0
                            if 0: print '@'
                            break
                        elif (begin > ds and end <= de) or\
                                (begin == ds and end < de):
                            dictrep[posset] += ':'+str(begin)+'^'
                            nooverlap = 0
                            if 0: print '^'
                            break
                        elif (begin <= ds and end >= de):
                            dictrep[posset] += ':'+str(begin)+'$'
                            nooverlap = 0
                            if 0: print '$'
                            break
                    #--------end tracing each position------------------
                    if nooverlap:
                        dictrep[posset] += ':'+str(begin)+'!'
                        if 0 : print '!'
                #----------if no domain-------------
                else:
                    dictrep[posset] += ':'+str(begin)+'!'
                    if 0 : print '!!'
                #--------------------------------
            #-----end of trace one group domain------------------
        #---------end of trace onelocusus-------------------------------
    #-------------end of all-------
    output(repDict)
#----------------end of function------------------
                        
if __name__ == '__main__':
    main()

