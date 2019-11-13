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
    This is designed to get the chemical structure from mass weight.

bg-file
ID  Formula Mass
CHEBI:90    C15H14O6    290.26810
CHEBI:598   C4H6O4R2    118.08800
CHEBI:776   C18H22O3    286.36550

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
from bisect import bisect_left
from re import findall
from copy import deepcopy

#from bs4 import BeautifulSoup

#reload(sys)
#sys.setdefaultencoding('utf8')

#from multiprocessing.dummy import Pool as ThreadPool

debug = 0

def fprint(content):
    print json_dumps(content,indent=1)

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
        metavar="FILEIN", help="A list of mass, one for each column")
    parser.add_option("-b", "--background-file", dest="bg",
        help="A three column file with format specified above.")
    parser.add_option("-D", "--diff", dest="diff",
        default=1.1, help="Maximum allowed differences for \
assign a formula to this mass. Default 1.1.")
    parser.add_option("-e", "--element", dest="element",
        default='', 
        help="Limit to given elements. Default all chemical compounds \
will be checked. Accept a list like \
<'C,H,O,N,P,S,K,Na,Br,F,Cl'>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readBg(bg, header=1):
    '''
    bgDict = {mass: [[ID, formula, mass], [ID2, formula2, mass2]]}
    '''
    bgDict = {}
    for line in open(bg):
        if header:
            header -= 1
            continue
        #--------------------------
        line = line.strip()
        id, formula, mass = line.split("\t")
        mass = float(mass)
        #assert mass not in bgDict, "Duplicate %f" % mass
        if mass not in bgDict:
            bgDict[mass] = [[id, formula, mass]]
        else:
            bgDict[mass].append([id, formula, mass])
    #-------------------------------
    bgMassL = bgDict.keys()
    bgMassL.sort()
    return bgMassL, bgDict
#-------------------------------
def findTheNearest(aim, searchSpace, start, end):
    '''
    This is designed to search for the closest number of <aim> in
    given <searchSpace>.
    
    Args:
        aim(number): number
        searchSpace(sorted list): numbers (searchX) to be searched
        start(number, 0-based): start poistion of searchSpace
        end(number, 0-based, not include): end position of searchSpace

    Returns:
        A set of the closest value and its position

    Raises:
        
    '''
    closest_pos = bisect_left(searchSpace, aim, start, end) - 1
    if closest_pos < 0:
        closest_pos = 0
    closest_value = searchSpace[closest_pos]
    diff = abs(aim - closest_value)
    if closest_pos+1 < end:
        diff2 = abs(aim - searchSpace[closest_pos+1])
        if diff > diff2:
            closest_pos = closest_pos + 1
            closest_value = searchSpace[closest_pos]
    #--------------------------------
    return closest_value,closest_pos
#---------------------------------------------------

def binarySearch(targetL, searchL, diff=0):
    '''
    This is designed to search the approximation number in <targetL>
    for each number in <searchL>.
    
    Args:
        targetL(sorted list): numbers (targetX) saved as targets
        searchL(sorted list): numbers (searchX) to be searched

    Returns:
        searchD(dict): {'searchX1':[targetX1, targetX2]}

    Raises:
        
    '''
    lo = 0
    hi = len(targetL)
    searchD = {}
    for num in searchL:
        num_left = num - diff
        num_right = num + diff
        pos_num_left  = bisect_left(targetL, num_left, lo, hi)
        lo = pos_num_left - 1
        if lo < 0: lo=0
        pos_num_right = bisect_left(targetL, num_right, lo, hi)+1
        if pos_num_right > hi:
            pos_num_right = hi

        tmpL = [(abs(targetL[i]-num), targetL[i]) \
            for i in range(lo, pos_num_right)]
        #print tmpL
        tmpL.sort(key=lambda x: x[0])
        searchD[num] = [i[1] for i in tmpL if i[0]<=diff]
    #-----------------------------------------------
    return searchD
#-----------------------------------------------
def parseFormula(formula):
    """ 
    Parse formula to element-level count

    Args:
        formula(string): like H2SO4

    Returns:
        formula(string): "HOS" (alphabet sorted)
        formulaD(dict): {'H':2, 'S':1, 'O':4}

    Raises:
    
    
    """
    formulaL = findall(r'([A-Z][a-z]*)(\d*)', formula)
    formulaD = []
    formulaS = set()
    for ele, cnt in formulaL:
        if cnt:
            cnt = int(cnt)
        else:
            cnt = 1
        formulaD[ele] = cnt
        formulaS.add(ele)
    #-------------------------------
    formulaL = list(formulaS)
    formulaL.sort()
    formula = ''.join(formulaL)
    return formula, formulaD
#-------------------------------------

def parseMassDformula(massD, elementD={}):
    for mass, matchLL in massD.items():
        removeL = []
        for mathL in matchL:
            formula = mathL[1]
            formula, formulaD = parseFormula(formula)
            delete = 0
            if elementD:
                for ele in formula:
                    if ele not in elementD:
                        delete = 1
                        break
            if delete:
                removeL.append(mathL)
            else:
                mathL.append(formula)
                mathL.append(formulaD)
        #----------------------
        for item in removeL:
            matchL.remove(item)
    #-----------------------------------
#--------------------------------
def sameElement(massDL_one, *massDL_other):
    '''
    This is designed to exclude formulas containing elements only in 
    <massDL_one>, but not in <massDL_other>.
    
    Args:
        massDL_one(2-dim list): [[ID, formula, mass, formula, formulaD], 
                         [ID2, formula2, mass2, formula2, formulaD2]]

        massDL_other(multiple 2-dim lists): 
                         [[ID, formula, mass, formula, formulaD], 
                         [ID2, formula2, mass2, formula2, formulaD2]]
    Returns:
        No return value

    Raises:
        
    '''
    local_massDL_one = deepcopy(massDL_one)
    local_one_formulaD = dict([(item_one, 1) \
            for item_one_L in local_massDL_one \
            for item_one in item_one_L[3]])

    local_massDL_otherL = deepcopy(massDL_other)
    local_other_formulaD = dict([(item_other, 1) \
            for local_massDL_other in local_massDL_otherL \
            for item_other_L in local_massDL_other \
            for item_other in item_other_L[3]])

    for item_one_L in local_massDL_one:
        for formula_ele in item_one_L[3]:
            if formula_ele not in local_other_formulaD:
                massDL_one.remove(item_one_L)
    #------------------------------------------------
    for local_massDL_other in local_massDL_otherL:
        for item_other_L in local_massDL_other:
            for item_other in item_other_L[3]:
                if item_other not in local_one_formulaD:
                    massDL_other.remove(item_other_L)
    #-------------------------------------------------------
#--------------------------------------


def integrateAnalysis(massL, massD, diff):
    '''
    This is designed to integrate two mass to screen the selection.
    The main principle is that if the sum of two masses equals to the
    other mass these three may have similar elements.

    We will first group mass and parse their chemical formula to get
    all elements of them.
    
    Args:
        massL(sorted list): numbers (massX) saved as targets
        massD(dict): {mass: [[ID, formula, mass, formula, formulaD], 
                         [ID2, formula2, mass2, formula2, formulaD2]]}

    Returns:
        conciseMassD(dict): {'searchX1':[targetX1, targetX2]}

    Raises:
        
    '''
    lenMassL = len(massL)
    for i in range(lenmassL-1):
        massI = massL[i]
        start = i
        for j in range(i+1, lenMassL):
            massJ = massL[j]
            totalM = massI + massJ
            clo_v, clo_p = findTheNearest(totalM, massI, start, lenMassL)
            start = clo_p - 1
            if abs(clo_v-totalM) <= diff:
                sameElement(massD[clo_v], massD[massI], massD[massJ])
            

#------integrateAnalysis-------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    bg   = options.bg
    diff = float(options.diff)
    massL = [float(i) for i in open(file)]
    massL.sort()
    elementD = dict([(i.strip(), 1) for i in options.element.split(',')])
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    bgMassL, bgDict = readBg(bg)
    searchD = binarySearch(bgMassL, massL, diff)
    massD = {}
    for mass in massL:
        targetL = searchD.get(mass, [])
        massD[mass] = []
        if targetL:
            for target in targetL:
                massD[mass].extend(bgDict[target])
            #--------------------------------------------
        else:
            #Here we will add a brute force crack
            print mass
    #-------------------------------------------------------
    parseMassDformula(massD, elementD)
    integrateAnalysis(massL, massD)



    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
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
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


