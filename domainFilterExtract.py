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

def extract(aDict, type, lenr=1):
    '''
    aDict = {locus: 
        {'overlap':[rep1, rep2],
        'new':[rep1, rep2],
        'dir':[rep1, rep2],
        'rid':[rep1, rep2]}
        }
    type : overlap, new, dir, rid
    lenr: a number
    newDict = {
        locus : [rep1, rep2]
    }
    '''
    assert isinstance(lenr, int)
    newDict = {}
    for locus, ldict in aDict.items():
        #locusOut = 0
        for key, value in ldict.items():
            if (key == type or (key != 'new' and type=='alloverlap')) \
                    and len(value):
                #if not locusOut:
                #    newDict[locus] = []
                #    locusOut = 1
                if locus not in newDict:
                    newDict[locus] = []
                for item in value:
                    if item.count(":") >= lenr:
                        newDict[locus].append(item)
                #----end rep output----------
            #-------end judging type---------
        #-------trace each locus
    #--------trace dict------------------------------
    keyL = newDict.keys()
    keyL.sort()
    for key in keyL:
        valueL = newDict[key]
        if len(valueL):
            print '>%s' % key
            print '\n'.join(newDict[key])
#--------------End extract---------------------------

def main():
    '''
    aDict = {locus: 
        {'overlap':[rep1, rep2],
        'new':[rep1, rep2],
        'dir':[rep1, rep2],
        'rid':[rep1, rep2]}
        }
    '''
    print >>sys.stderr, "Print the result to screen"
    lenargv = len(sys.argv)
    if lenargv != 2 and lenargv != 4:
        print >>sys.stderr, 'Using python %s filename [choice<new,\
overlap, dir, rid, alloverlap>] [reptime<1,2,3..>]. If only the first parameter\
given, it will output al results.'\
            % sys.argv[0]
        sys.exit(0)
    #-----------------------------------------
    aDict = {}
    for line in open(sys.argv[1]):
        if line[0] == '>':
            locus = line[1:-1]
            aDict[locus] = {}
            (aDict[locus])['overlap'] = [] #@
            (aDict[locus])['new'] = [] #!
            (aDict[locus])['dir'] = [] #$
            (aDict[locus])['rid'] = [] #^
        else:
            assert line.count('!')+line.count('@')+line.count('$')\
                +line.count('^') == line.count('#')+1
            assert line.count('#')+1 == line.count(':')
            lineL = line[:-1].split('#')
            overlap = ''
            new = ''
            dir = ''
            rid = ''
            for item in lineL:
                if item[-1] == '!':
                    new += item.replace('!', '#')
                elif item[-1] == '@':
                    overlap += item.replace('@', '#')
                elif item[-1] == '^':
                    rid += item.replace('^', '#')
                elif item[-1] == '$':
                    dir += item.replace('$', '#')
                else:
                    print 'Wrong end'
                    sys.exit(1)
            #-------after trace------------------
            if len(new) > 0:
                (aDict[locus])['new'].append(new[:-1]) #!
            if len(overlap) > 0:
                (aDict[locus])['overlap'].append(overlap[:-1]) #@
            if len(dir) > 0:
                (aDict[locus])['dir'].append(dir[:-1]) #$
            if len(rid) > 0:
                (aDict[locus])['rid'].append(rid[:-1]) #^
        #-----------------finish else------------------
    #----------finish reading-------------
    if lenargv == 2:
        #--------save the number of different protein and repetition
        pDict = {'overlap':0, 'new':0, 'dir':0, 'rid':0}
        rDict = {'overlap':0, 'new':0, 'dir':0, 'rid':0}
        #--------this demands rep times larger than 3
        pDictL = {'overlap':0, 'new':0, 'dir':0, 'rid':0}
        rDictL = {'overlap':0, 'new':0, 'dir':0, 'rid':0}
        '''
        aDict = {locus: 
            {'overlap':[rep1, rep2],
            'new':[rep1, rep2],
            'dir':[rep1, rep2],
            'rid':[rep1, rep2]}
            }
        '''
        for locus, ldict in aDict.items():
            for key, value in ldict.items():
                repN = len(value)
                if repN:
                    pDict[key] += 1
                    rDict[key] += repN
                    py = 0 # a flag to see if this protein is ok
                    for item in value:
                        if item.count(':') > 2:
                            rDictL[key] += 1
                            py = 1
                    if py:
                        pDictL[key] += 1
        #----------------------------------------------
        filel= [sys.argv[1]+'.pn', sys.argv[1]+'.rn',\
                sys.argv[1]+'.pn3', sys.argv[1]+'.rn3']
        dictl= [pDict, rDict, pDictL, rDictL]
        keyL = ['new','overlap','dir', 'rid']
        for i in range(4):
            fh = open(filel[i], 'w')
            for key in keyL:
                print >>fh, '\t'.join((key, str((dictl[i])[key])))
            fh.close()
    elif lenargv == 4:
        extract(aDict, sys.argv[2], int(sys.argv[3]))
if __name__ == '__main__':
    main()

