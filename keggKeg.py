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
import re

def func_patb(line, patb0, patb1, patb2, patb3, patb4,
    patb, secDict, aDict):
    matchb = patb0.match(line)
    if matchb:
        b = matchb.group(1)
        b = b.strip()
        for tmpDict in aDict:
            if isinstance(tmpDict, dict):
                if b in tmpDict:
                    print >>sys.stderr, "Wrong, duplicate B", b
                    sys.exit(1)
        #----end for-------------------------------
        #b = b.replace(' ', '_')
        #b = b.replace('/', '_')
        secDict[b] = []
        aDict.append(secDict)
        return b
    else:
        matchb = patb1.match(line)
        if matchb:
            b = '\t'.join((matchb.group(1), matchb.group(2),
                matchb.group(3), matchb.group(5), matchb.group(4)))
            #assert b not in aDict, b
            aDict.append(b)
            return 0
        else:
            matchb = patb2.match(line)
            if matchb:
                b = '\t'.join((matchb.group(1), matchb.group(2),
                    matchb.group(3), 'EC:', matchb.group(4)))
                #assert b not in aDict, b
                aDict.append(b)
                return 0
            else:
                matchb = patb3.match(line)
                if matchb:
                    b = '\t'.join(('NULL', matchb.group(1),
                        matchb.group(2), 'EC:', matchb.group(3)))
                    #assert b not in aDict, b
                    aDict.append(b)
                    return 0
                else:
                    matchb = patb4.match(line)
                    if not matchb:
                        matchb = patb.match(line)
                    b = matchb.group(1)
                    b = b.strip()
                    for tmpDict in aDict:
                        if isinstance(tmpDict, dict):
                            if b in tmpDict:
                                print >>sys.stderr, "Wrong, duplicate B", b
                                sys.exit(1)
                    #----end for-------------------------------
                    #b = b.replace(' ', '_')
                    #b = b.replace('/', '_')
                    secDict[b] = []
                    aDict.append(secDict)
                    return b
                #--------------End non-regular class B-------- 
            #---------ENd no locus-----label by NULL-----
        #----------End no-enzyme-----------
    #----------End enzyme--------------
#---------End func_patb-------------
def readKeg(file, aDict):
    '''
    aDict = {A:[{B:[{C:[locus\tKO\tDES\tEC\tEName, ]}]}]}
    A: pathway based gene classification
    B: protein family classification
    C: Compound classification
    D: drug classification
    '''
    #Before I think B-lines have no item only a definition, but it is
    #not true. 
    patb = re.compile(r"B {2}(.+)")
    patb0 = re.compile(r"B {2}<b>(.*?)</b>")
    patb1 = re.compile(r"B {2}([^;]+); <.+?>(.+)<.+?> +([^;]+); (.+?) \[(.+)\]") #3
    patb2 = re.compile(r"B {2}([^;]+); <.+?>(.+)<.+?> +([^;]+); (.+)") #4
    patb3 = re.compile(r"B {2}<.+?>(.+)<.+?> +([^;]+); (.+)") #5
    patb4 = re.compile(r"B {2}<.+?>.*?<.+?> +(.+)") #1 # fifth
    #-----changed 2011-12-21-------------------
    #Before I used #1 to prase file q00001.key. Most 'compound
    #classification' which have offsprings all followed by '[PATH]'.
    #Later on, I find a suffix '[BR] by chance'. So I changed 'patc'
    #to #2.
    #patc = re.compile(r"C {4}(\d+) (.+?) \[PATH") #1
    #Also, there may have gene locus here as in D.
    patc = re.compile(r"C {4}([^[]+)") #2 
    patc0 = re.compile(r"C {4}<b>(.*?)</b>")
    patc1 = re.compile(r"C {4}([^;]+); <.+?>(.+)<.+?> +([^;]+); (.+?) \[(.+)\]") #3
    patc2 = re.compile(r"C {4}([^;]+); <.+?>(.+)<.+?> +([^;]+); (.+)") #4
    patc3 = re.compile(r"C {4}<.+?>(.+)<.+?> +([^;]+); (.+)") #5
    patc4 = re.compile(r"C {4}<.+?>.*?<.+?> +(.+)") #1 # fifth
    #-------recorded 2011-12-21--------------------
    #In #1, patd is the complete form, patd2 is used to deal with
    #lines without enzyme nominate. 
    #Unexpectedly, there is a third type of D-starting lines, which
    #have only one ';' and no gene locus.
    patd = re.compile(r"D {6}([^]]+)") #2 # sixth
    patd0 = re.compile(r"D {6}<b>(.*?)</b>") #first
    patd1 = re.compile(r"D {6}([^;]+); <.+?>(.+)<.+?> +([^;]+); (.+?)\[(.+)\]") #1 #second
    patd2 = re.compile(r"D {6}([^;]+); <.+?>(.+)<.+?> +([^;]+); (.+)") #1 #third
    patd3 = re.compile(r"D {6}<.+?>(.+)<.+?> +([^;]+); (.+)") #1 #forth
    patd4 = re.compile(r"D {6}<.+?>.*?<.+?> +(.+)") #1 # fifth
    #---------add 2011-12-21-------------------------
    pate = re.compile(r"E {8}([^]]+)") #2 
    pate0 = re.compile(r"E {8}<b>(.*?)</b>")
    pate1 = re.compile(r"E {8}([^;]+); <.+?>(.+)<.+?> +([^;]+); (.+?) \[(.+)\]") #1
    pate2 = re.compile(r"E {8}([^;]+); <.+?>(.+)<.+?> +([^;]+); (.+)") #1
    pate3 = re.compile(r"E {8}<.+?>(.+)<.+?> +([^;]+); (.+)") #1
    pate4 = re.compile(r"E {8}<.+?>.*?<.+?> +(.+)") #1 # fifth
    #---------add 2011-12-21-------------------------
    patf = re.compile(r"F {10}([^]]+)") #2 
    patf0 = re.compile(r"F {10}<b>(.*?)</b>")
    patf1 = re.compile(r"F {10}([^;]+); <.+?>(.+)<.+?> +([^;]+); (.+?) \[(.+)\]") #1
    patf2 = re.compile(r"F {10}([^;]+); <.+?>(.+)<.+?> +([^;]+); (.+)") #1
    patf3 = re.compile(r"F {10}<.+?>(.+)<.+?> +([^;]+); (.+)") #1
    patf4 = re.compile(r"F {10}<.+?>.*?<.+?> +(.+)") #1 # fifth

    for line in open(file):
        if line.startswith('#ENTRY'):
            ko = line.strip().split()[1]
        elif line.startswith('#DEFINITION'):
            ko += '_'+'_'.join(line.split()[1:-1])
        elif line.startswith('A'):
            a = line[4:-5] 
            if a in aDict:
                print >>sys.stderr, "Wrong, duplicate A"
            aDict[a] = []
        elif line.startswith('B  '):
            secDict = {}
            b = func_patb(line, patb0, patb1, patb2, patb3, 
                    patb4, patb, secDict, aDict[a])
        elif line.startswith("C    "):
            assert b != 0
            thirdDict = {}
            c = func_patb(line, patc0, patc1, patc2, patc3, patc4,
                    patc, thirdDict, secDict[b])
        elif line.startswith("D      "):
            assert c != 0
            forthDict = {}
            d = func_patb(line, patd0, patd1, patd2, patd3, patd4
                    ,patd, forthDict, thirdDict[c])
        elif line.startswith("E        "):
            assert d != 0
            fifthDict = {}
            e = func_patb(line, pate0, pate1, pate2, pate3, pate4,
                    pate, fifthDict, forthDict[d])
        elif line.startswith("F          "):
            assert e != 0
            sixthDict = {}
            f = func_patb(line, patf0, patf1, patf2, patf3, patf4,
                    patf, sixthDict, fifthDict[e])
            assert f == 0
        else:
            print >>sys.stderr, line,
    #-----------------------------------------------------------
    return ko
#-----------------------------------------------------------

def output(aDict, ko, fhp):
    '''
    aDict = {A:[{B:[{C:[locus\tEC\tKO\tDES, ]}]}]}
    A: pathway based gene classification
    B: protein family classification
    C: Compound classification
    D: drug classification
    aDict is a dict, it use A as key, its value is a list. The element
    of the list may be a str or a dict like aDict.
    ko is the head of the file recording the ENTRY+DEFINITION.
    #*fhp is a tuple of file handles.
    fhp is a list of file handlers.
    '''
    fhDict = {}
    aKeyL = aDict.keys()
    aKeyL.sort()
    for a in aKeyL:
        prefix = (ko+'.'+a).replace(' ', '_')
        prefix = prefix.replace('/', '_')
        prefix = prefix.replace('(', '_')
        prefix = prefix.replace(')', '_')
        prefix = prefix.replace('___', '_')
        prefix = prefix.replace('__', '_')
        prefix = prefix.strip('_')
        #------------------------------
        valueB = aDict[a]
        if len(valueB) == 0:
            print >>sys.stderr, prefix
            continue
        #------------------------------------
        fhDict[a] = fhp[0:]
        fh = open(prefix, 'w')
        fhDict[a].append(fh)
        #--------------------------------
        #Here, I assert everything in the list have same attribute,
        #which means they are either all dict or all string.
        if isinstance(valueB[0], dict):
            for potentialBDict in valueB:
                output(potentialBDict, prefix, fhDict[a][0:])
        elif isinstance(valueB[0], str):
            for fh2 in fhDict[a]:
                try:
                    print >>fh2, '\n'.join(valueB)
                except:
                    print '----------'+str(fh2)+'--------'
                    print '\n'.join(valueB)
                    print '----------'+str(fh2)+'--------'
        fh.close()
#-------------------------------------------------

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) < 2:
        print >>sys.stderr, 'Using python %s filename' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------------
    aDict = {}
    '''
    aDict = {A:[{B:[{C:[locus\tEC\tKO\tDES, ]}]}]}
    A: pathway based gene classification
    B: protein family classification
    C: Compound classification
    D: drug classification
    '''
    #for file in sys.argv:
    file = sys.argv[1]
    ko=readKeg(file, aDict)
    output(aDict, ko, [])
if __name__ == '__main__':
    main()

