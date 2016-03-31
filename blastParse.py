#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
from optparse import OptionParser as OP
import decimal
import sys
import re
import os

DEBUG       = 0
TYPE        = 1 #means simple
IDENTITY    = 1 #means identity
UNIQ        = 0 #locus or segments
ALL         = 0 #output all info 

def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python '+ argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    desc="Prase Blast Output\n\n"
    usages = desc+"%prog [-i blast] [-o output]"
    parser = OP(usage=usages)
    parser.add_option("-f", "--input-file", dest="filein",
        metavar="FILEIN", help="The result of blast default output")
    parser.add_option("-o", "--output-file", dest="fileout",
        metavar="FILEOUT", help="Save the parse [stdout]")
    parser.add_option("-p", "--positive", dest="chose", default=True,
        help="Using the positive value [default identity]", 
        action="store_false")
    parser.add_option("-v", "--value", dest="value",
        nargs = 2, default=(0.7, 1.0), type="float",
        help="The min and max value seted [%default]")
    parser.add_option("-c", "--complex", dest="type", default=True,
        help="Complex mode [Simple, no merge]", action="store_false")
    parser.add_option("-d", "--debug", dest="verbose", default=False,
        help="Open the debug mode [False]", action="store_true")
    parser.add_option("-u", "--uniq", dest="uniq", default=False,
        help="Output each locus or each segments [all segments]",
        action="store_true")
    parser.add_option("-a", "--all", dest="out", default=False,
        action="store_true", 
        help="Output all the informatin, include query locus, querydescription,\
query length, hit locus, hit description,hit-locus-length,matchLen,\
matchPosition [no desc, no pos]")
    #parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv[1:])
    if options.filein == None:
        print >>sys.stderr, "Wrong filename after -f"
        sys.exit(1)
    return (options, args)
#--------------------------------------------------------------------

def extractInfoFromResult(filein, aDict):
    '''
    aDict = 
        {(queryLoc, queryLen)
        :
            {(hitLoc, hitLen)
            :
                {(identity,matchLen,positive,label)
                :
                [(queryPos),(sbjctPos)]//before tune
                (queryPos1, queryPos2, sbjctPos1, sbjctPos2)//after tune
                }
            }
        }
    '''
    global DEBUG
    debug = DEBUG
    queryLoc = re.compile('Query= (.+)')
    queryLen = re.compile('.*?(\d+) letters')
    hitLoc   = re.compile('>.+?\|+(.+?)[\|\s]')
    hitLen   = re.compile('\s+Length = (\d+)') 
    identity = re.compile(' Identities = (\d+?)\/(\d+).*?Positives = (\d+?)\/')
    queryPos = re.compile('Query: (\d+?)\s*?.+?\s+(\d+)')
    sbjctPos = re.compile('Sbjct: (\d+?)\s*?.+?\s+(\d+)')
    #---------Begin Reading------------------------------------------
    #---------------some flags---------------------------------------
    fh = open(filein)
    while 1:
        lamb = 0
        #-------------Get queryLoc-----------------------------------
        while 1:
            line = fh.readline()
            if debug: print >>sys.stderr, 'queryLoc while: %s*' % line[:-1]
            if len(line) == 0: break
            if line.find('Lambda') != -1: 
                lamb = 1
                break
            if line[:6] == 'Query=':
                queryLocV = queryLoc.match(line).group(1)
                if debug: print >>sys.stderr, queryLocV
                break
        #-------------Get queryLen-----------------------------------
        if len(line) == 0 or lamb: break
        while 1:
            line = fh.readline()
            assert(len(line) != 0)
            if line.find('letters') != -1:
                break
        queryLenV = int(queryLen.match(line).group(1))
        if debug: print >>sys.stderr, queryLenV
        keySet = (queryLocV, queryLenV)
        aDict[keySet] = {}
        #---------------if hitted--------------------------------
        nohit = 0
        while 1:
            line = fh.readline()
            if debug: print >>sys.stderr, 'if hit while: %s*' % line[:-1]
            assert(len(line) != 0)#: break
            if line.find("No hits found") != -1:
                nohit = 1
                break
            elif line.find("Sequences producing significant alignments") != -1:
                break
        if len(line) == 0: break
        if nohit: continue
        #-------------Get hit---------------------------------------
        hitFlag = 0
        while 1:
            #-------------Get hitLoc-------------------------------------
            while 1:
                if debug: print >>sys.stderr, 'hitLoc while: %s*' % line[:-1]
                if line[0] != '>':
                    line = fh.readline()
                #if len(line) == 0: break
                assert(len(line) != 0)
                if line[0] == '>':
                    hitLocV = hitLoc.match(line).group(1)
                    if debug:
                        print >>sys.stderr, 'hitLoc while: %s*' % line[:-1]
                        print >>sys.stderr, 'hitLocV', hitLocV
                    break
            if len(line) == 0: break
            #------------Get hitLen--------------------------------------
            while 1:
                line = fh.readline()
                if debug: print >>sys.stderr, 'hitLen while: %s*' % line[:-1]
                #if len(line) == 0: break
                assert(len(line) != 0)
                if line.find("Length = ") != -1:
                    hitLenV = int(hitLen.match(line).group(1))
                    if debug: print >>sys.stderr, 'hitLenV', hitLenV
                    break
            hitKeySet = (hitLocV, hitLenV)
            aDict[keySet][hitKeySet] = {}
            #------------Get identity, query pos, sbjct pos--
            while 1:
                #------------Get each identity, query pos, sbjct pos--
                label = 0
                while 1:
                    line = fh.readline()
                    if debug: print >>sys.stderr, 'identity while: %s*' % line[:-1]
                    #if len(line) == 0: break
                    assert(len(line) != 0)
                    if line.find('Identities') != -1:
                        identityV, matchLen, positive = \
                            identity.match(line).group(1,2,3)
                        identityKeySet = (int(identityV), int(matchLen), 
                            int(positive), label)
                        label += 1
                        if debug: print >>sys.stderr, identityKeySet
                        aDict[keySet][hitKeySet][identityKeySet] = []
                        break
                if len(line) == 0: break
                while 1:
                    line = fh.readline()
                    length = len(line)
                    assert(length != 0)
                    if length == 1: break
                while 1:
                    line = fh.readline()
                    if debug: print >>sys.stderr, 'query while: %s*' % line[:-1]
                    assert(len(line) != 0)
                    if line == '\n': break
                    queryPosV1, queryPosV2 = \
                        queryPos.match(line).group(1,2)
                    fh.readline() #throw the mid line
                    sbjctPosV1, sbjctPosV2 = \
                        sbjctPos.match(fh.readline()).group(1,2)
                    if debug:
                        print >>sys.stderr, queryPosV1, queryPosV2, sbjctPosV1, sbjctPosV2 
                    fh.readline() #throw the blank line
                    aDict[keySet][hitKeySet][identityKeySet].\
                        append((int(queryPosV1), int(queryPosV2)))
                    aDict[keySet][hitKeySet][identityKeySet].\
                        append((int(sbjctPosV1), int(sbjctPosV2)))
                line = fh.readline()
                if debug: print >>sys.stderr, '177: %s*' % line[:-1]
                if line == '\n': continue
                else:
                    if line[0] != '>': hitFlag = 1
                    break
                #----------End Get each--------------------
            #--------------End identity, query pos, sbjct pos--------
            if hitFlag: break
        #------------------End of get hit----------------------------
    #----------------------End of reading----------------------------
    fh.close()
    if debug: print >>sys.stderr, aDict
    
    tuneAdict(aDict)
    if debug:
        print >>sys.stderr, aDict
#--------------------------------------------------------------------

def tuneAdict(aDict):
    '''
    aDict = 
        {(queryLoc, queryLen)
        :
            {(hitLoc, hitLen)
            :
                {(identity,matchLen,positive, label)
                :
                [(queryPos),(sbjctPos)]//before tune
                (queryPos1, queryPos2, sbjctPos1, sbjctPos2)//after tune
                }
            }
        }
    '''
    for hitDict in aDict.values():
        for idenDict in hitDict.values():
            for qbkey, qbval in idenDict.items():
                idenDict[qbkey] = (qbval[0][0], qbval[-2][1], \
                    qbval[1][0], qbval[-1][1])
#-----------End of tuneAdict-----------------------------------------
def writeDict(aDict, fileout):
    '''
    '''
    if fileout != None:
        fh = open(fileout, 'w')
        oldstd = sys.stdout
        sys.stdout = fh
    for key, hitDict in aDict.items():
        print key
        for hitlocus, idenDict in hitDict.items():
            print '\t', hitlocus
            for qbkey, qbval in idenDict.items():
                print '\t\t', qbkey
                print '\t\t\t', qbval
    if fileout != None:
        fh.close()
        sys.stdout = oldstd
#----------End of writeDict------------------------------------------    
def merge(aDict):
    pass
#----------End of merge----------------------------------------------

def getwanted(aDict, minV, maxV, fileout):
    '''
    aDict = 
        {(queryLoc, queryLen)
        :
            {(hitLoc, hitLen)
            :
                {(identity,matchLen,positive,label)
                :
                [(queryPos),(sbjctPos)]//before tune
                (queryPos1, queryPos2, sbjctPos1, sbjctPos2)//after tune
                }
            }
        }
    '''
    ctDecimal = decimal.Decimal
    global IDENTITY
    global UNIQ
    global ALL
    uniq = UNIQ
    all = ALL
    if IDENTITY: #means using identity value
        index = 0
    else:
        index = 2
    resultDict = {}
    for key, hitDict in aDict.items():
        resultDict[key] = []
        length = key[1]
        for hitlocus, idenDict in hitDict.items():
            idenDictKey = idenDict.keys()[0]
            identity = idenDictKey[index]
            idenPercent = identity / length
            idenPercent = ctDecimal(str(idenPercent))
            if minV <= idenPercent <= maxV:
                if all:
                    newse = (hitlocus, identity, idenDict[idenDictKey])
                else:
                    newse = (hitlocus, identity)
                resultDict[key].append(newse)
                if uniq: break
    #-----------END for----------------------------------------------
    if fileout != None:
        fh = open(fileout, 'w')
        oldstd = sys.stdout
        sys.stdout = fh
    #print resultDict
    for key, value in resultDict.items():
        if len(value) == 0: continue
        print '>%s\t%d' % (key[0], key[1])
        for item in value:
            if all:
                print '\t%s\t%d\t%d\tquery:(%d, %d),sbjct:(%d, %d)' %\
                    (item[0][0],item[0][1], item[1], \
                    item[2][0], item[2][1], item[2][2], item[2][3])
            else:
                print '\t%s\t%d\t%d' % (item[0][0], item[0][1], item[1])
    if fileout != None:
        fh.close()
        sys.stdout = oldstd
#----------End of getwanted------------------------------------------



def main():
    (options, args) = cmdparameter(sys.argv)
    #----------------------------------------------------------------
    filein = options.filein
    fileout = options.fileout
    global DEBUG
    global TYPE
    global IDENTITY
    global UNIQ
    global ALL
    DEBUG = options.verbose
    TYPE = options.type #True means simple, False means complex
    IDENTITY = options.chose #True means identity, false means positive
    UNIQ = options.uniq
    ALL = options.out
    minV, maxV = options.value
    minV = decimal.Decimal(str(minV))
    maxV = decimal.Decimal(str(maxV))
    aDict = {}
    extractInfoFromResult(filein, aDict)
    if DEBUG: writeDict(aDict, fileout)
    if not TYPE:
        merge(aDict)
    getwanted(aDict, minV, maxV, fileout)
    


if __name__ == '__main__':
    main()

