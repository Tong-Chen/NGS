#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to simulate Mass Spectrum in silicon.

Input file format given to `-f`

****************************************************************
# Lines begin with pound sign (#) are used for comments and will
# be ignored when reading in programs.
# Blank lines will be ignored too.
# Below lists the format of input
# 取代基排列模式: 可以为`sequential`或`permutation`
# sequential: 表示取代基按给定的顺序连接到母核上
# permutation: 表示取代基之间先做个排列组合再连到母核上
# 母核
# <tab>取代基1:取代基最大数目<tab>取代基排列模式<tab>电荷模式
# <tab>取代基1:取代基最大数目;取代基2:取代基最大数目;<tab>取代基排列模式<tab>电荷模式

#Compound 1
C6
	H2O:3   sequential  -
	C5:3;O5:3   permutations    +

#Compound 2
C5	
	C5:2;O5:2   sequential  +
	H2O:3   permutations    -

****************************************************************

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from chempy import Substance
from_formula = Substance.from_formula


#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')
from itertools import permutations

debug = 0

def readInput(input):
    """
    inputD = {'mother_nucleus': [
                [[(R1, 3)], 'sequential', '+'], 
                [[(R2, 3), (R3, 3)], 'permutation', '-']
               ], 
             }
    """
    inputD = {}
    for line in open(input):
        if line.startswith('#') or line.strip() == "":
            continue
        if line.startswith('\t'):
            line, type = line.strip().split('\t')
            lineL = line.split(';')
            subst_grpL = [item.split(':') for item in lineL]
            subst_grpL = [subst_grpL, type]
            inputD[key].append(subst_grpL)
        else:
            key = line.strip()
            assert key not in inputD, "Duplicate key %s" % key
            inputD[key]= []
    return inputD
#--------------------------------------

def combine(inputD):
    """
    inputD = {'mother_nucleus': [
                [[(R1, 3)], 'sequential'], 
                [[(R2, 3), (R3, 3)], 'permutation']
               ], 
             }
    subgrpLL = [ [[(R1,3)],'sequential'], 
                 [[(R2,3),(R3,3)], 'permutation'] ]

    subgrpL  = [[ (R2,3), (R3,3) ], 'permutation']

    compoundL = [[mother_nucleus,R1],
                    [mn,R1,R1],[mn,R1,R1,R1], 
                 [mn,R2,R3], [mn,R2,R3,R2,R3], 
                     [mn,R2,R3,R2,R3,R2,R3]
                ]

    """
    compoundL = []
    for mn, subgrpLL in inputD.items():
        if debug:
            print >>sys.stderr, "##Mother nuclear\tsubgrpLL"
            print >>sys.stderr, mn, subgrpLL
        for subgrpL in subgrpLL:
            type = subgrpL[1]
            grpL = [i[0] for i in subgrpL[0]]
            if type == "permutation": 
                combGrpL = list(permutations(grpL, len(grpL)))
            elif type == "sequential":
                combGrpL = [grpL]
            else:
                print >>sys.stderr, "Unsupported type: %s" % type
                sys.exit(1)
            count = int(subgrpL[0][0][1])
            if debug:
                print >>sys.stderr, "##combGrpL"
                print >>sys.stderr, combGrpL
            for combGrp in combGrpL:
                #combGrp = list(combGrp) * count
                for cnt in range(count):
                    cnt += 1
                    tmpL = [mn]
                    tmpL.extend(list(combGrp)*cnt)
                    compoundL.append(tmpL[:])
                if debug:
                    print >>sys.stderr, "##compoundL"
                    print >>sys.stderr, compoundL
    return compoundL
#--------------------------------------

def massCompute(formulaAddL, formulaSubstractL=[]):
    """
    Use chempy.Substance.from_formula to get molecular mass of given
    chemical formulas.

    Args:
        formulaAddL(list): The mass of chemical formulas in this param
                            will be added together.
        formulaSubstractL(list): The mass of chemical formulas in this param
                            will be substracted from above summarized
                            mass.
    """
    mass = 0
    for formula in formulaAddL:
        cpd = from_formula(formula)
        mass += cpd.mass

    for formula in formulaSubstractL:
        cpd = from_formula(formula)
        mass -= cpd.mass

    return mass
#------massCompute----------------------------


def massWhole(compound, massSubCpdL):
    """
    This is used to get the mass of whole molecular.

    Args:
        ##compound(list):[chemial1, chemical2, ...]
        compound(list):   [mother_nucleus,R1],
                       or [mn,R1,R1]
                       or [mn,R2,R3] 
                       or [mn,R2,R3,R2,R3]

    Process:
        1. Summarized mass of all parts and minus H2O
        2. Result of 1 + COOH - H2O
        3. Result of 1 * 2 - 1
        4. Result of 1 + Cl
    """
    cpd = ';'.join(compound)
    len_cpd = len(compound)
    mass_cpd = massCompute(compound, ["H2O"]*(len_cpd-1))
    massSubCpdL.append([cpd, mass_cpd])

    ###Add cooh
    cpd_cooh = cpd+';'+'cooh'
    mass_cpd_cooh = mass_cpd + massCompute(["COOH"])
    massSubCpdL.append([cpd_cooh, mass_cpd_cooh])

    ###Add Cl
    cpd_Cl = cpd+';'+'Cl'
    mass_cpd_Cl = mass_cpd + massCompute(["Cl"])
    massSubCpdL.append([cpd_Cl, mass_cpd_Cl])

    ###2 * cpd - 1
    cpd2 = cpd + " * 2"
    mass_cpd2 = 2 * mass_cpd - 1
    massSubCpdL.append([cpd2, mass_cpd2])


#---massWhole---------------------------------

def massSimulate(compound, massCpdD):
    """
    This is used to simulate the MS process by splitting compounds.

    Args:
        compound(list):   [mother_nucleus,R1],
                       or [mn,R1,R1]
                       or [mn,R2,R3] 
                       or [mn,R2,R3,R2,R3]

        massCpdD(dict): {cpd: [
                            [part1, mass1], 
                            [part2, mass2]
                         ]}                    

    """ 
    #Compute the mass of whole molecular
    cpd_key = ';'.join(compound)
    assert cpd_key not in massCpdD, "Duplicate %s" % cpd
    massSubCpdL = []
    massWhole(compound, massSubCpdL)
    
    #Get other combinations
    cpdL = []
    cntCpd = len(compound)
    for i in range(cntCpd):
        if i:
            end = cntCpd + 1
        else:
            end = cntCpd
        for j in range(i+1, end):
            tmpCpd = compound[i:j]
            if tmpCpd not in cpdL:
                cpdL.append(tmpCpd)
    
    #Compute the mass for each part
    for cpd in cpdL:
        cpdStr = ';'.join(cpd)
        len_cpd = len(cpd)
        substractH = ['H']
        substractH.extend(["H2O"]*(len_cpd-1))
        mass = massCompute(cpd, substractH)
        if debug:
            print >>sys.stderr, massSubCpdL
            print >>sys.stderr, cpdStr
            print >>sys.stderr, mass
        massSubCpdL.append([cpdStr+';-H', mass])
        substractH2O = ["H", "H2O"]
        substractH2O.extend(["H2O"]*(len_cpd-1))
        mass = massCompute(cpd, substractH2O)
        massSubCpdL.append([cpdStr+';-H-H2O', mass])
    #--------------------------------------------------
    massCpdD[cpd_key] = massSubCpdL[:]

#----massSimulate----------------------------------

def massSplit(compoundL, aDict):
    """
    This is used to simulate the MS process by splitting compounds.

    Args:
        compoundL(list):   [mother_nucleus,R1],
                       or [mn,R1,R1]
                       or [mn,R2,R3] 
                       or [mn,R2,R3,R2,R3]

    Desp:
        Take compoundL <[mn,R2,R3]> as an example
        It will be split into:
            [mn], [R2,R3]
            [mn], [R2], [R3]
            [mn, R2], [R3]
        #------------------------------------------
        Take compoundL <[mn,R2,R3,R4]> as an example
        It will be split into:
            [mn], [R2,R3,R4]
            [mn], [R2], [R3, R4]
            [mn], [R2], [R3], [R4]
            [mn], [R2, R3], [R4]
            [mn, R2], [R3, R4]
            [mn, R2], [R3], [R4]
            [mn, R2, R3], [R4]
            
    """ 
    len_compoundL = len(compoundL) - 1
    for i in range(len_compoundL):
        firstL  = compoundL[0:i+1]
        key = ';'.join(firstL)
        secondL = compoundL[i+1:]
        value = ';'.join(secondL)
        aDict[key] = {}
        aDict[key][value] = i
        if len(secondL) > 1:
            massSplit(secondL, aDict[key])
#-------END massSplit------------------------------

def transferSplitDict(aDict, cpdL, tmpL):
    """
    input: aDict = 
        {mn: 
            {R2: 
                {
                R3: {R4: 0},  
                R3;R4: 0
                }, 
             R2;R3: 
                {R4: 1}, 
             R2;R3;R4: 0}, 

          mn;R2: 
             {R3: 
                {R4: 0},  
                 R3;R4: 1}, 

          mn;R2;R3: 
                {R4: 2}
        }

    output: cpdL =
        [[mn,R2,R3,R4],[mn,R2,R3;R4],[mn,R2;R3,R4],[mn,R2;R3;R4],
         [mn;R2,R3,R4],[mn;R2,R3;R4],
         [mn;R2;R3,R4]
        ]
    
    tmpL: a temporary list used for iterating
    """
    if debug:
        print >>sys.stderr, "Input aDict:", aDict
        print >>sys.stderr, "Input tmpL", tmpL
        print >>sys.stderr, "Input tmpL id", id(tmpL)
        print >>sys.stderr, "Input cpdL", cpdL

    for key, value in aDict.items():
        if debug:
            print >>sys.stderr, "++Current key", key
        ## Introducing innertmpL is the main point.
        ## It keeps the original version of tmpL which should be 
        ## shared by all items in given dict.
        innertmpL = tmpL[:]
        innertmpL.append(key)
        if type(value) != dict:
            cpdL.append(innertmpL[:])
            if debug:
                print >>sys.stderr, "\tReturnned value", cpdL
        else:
            if debug:
                print >>sys.stderr, "\tdict for transferSplitDict", value
                print >>sys.stderr, "\ttmpL for transferSplitDict", tmpL
                print >>sys.stderr, "\ttmpL for transferSplitDict id", id(tmpL)
                print >>sys.stderr, "\tcpdL for transferSplitDict", cpdL
            transferSplitDict(value, cpdL, innertmpL[:])
    #------------END for-----------------------
#-------------transferSplitDict----------------------------------

def massSimulate2(compoundL, massCpdD):
    """
    This is used to simulate the MS process by splitting compounds.

    Args:
        compoundL(list):   [mother_nucleus,R1],
                       or [mn,R1,R1]
                       or [mn,R2,R3] 
                       or [mn,R2,R3,R2,R3]

        massCpdD(dict): {cpd: {
                              whole: [(mode1, mass1), (mode2, mass2), ]
                              cpd_split1: [(part1_1,mass1),(part1_2,mass2),],
                              cpd_split2: [(part2_1,mass1),(part2_2,mass2),],
                            }, 
                        } 
    """ 
    #Used to record already computed mass for each compound only.
    #Put this variable at the beginning of the script will make
    #it as a global variable.
    massRecord = {}
    #Compute the mass of whole molecular
    cpd_key = ';'.join(compoundL)
    assert cpd_key not in massCpdD, "Duplicate %s" % cpd
    massSubCpdL = []
    massWhole(compoundL, massSubCpdL)
    massCpdD[cpd_key] = {}
    massCpdD[cpd_key]['whole'] = massSubCpdL
    
    #Get other combinations
    cpdL = []
    compoundDict = {}
    massSplit(compoundL, compoundDict)
    transferSplitDict(compoundDict, cpdL, [])
    #sort unneeded, since all variables will be save in a dict
    #do the sort when output
    #cpdL.sort()

    #Compute the mass for each part
    #cpdL = [  [mn,R2,R3,R4],[mn,R2,R3;R4],
    #          [mn,R2;R3,R4],[mn,R2;R3;R4]
    #       ]

    '''
    # Every cpd can be split in multiple ways and would be treated
    # as a list.

    massCpdD(dict): {cpd: {
                          whole: [(mode1, mass1), (mode2, mass2), ]
                          cpd_split1: [(part1_1,mass1),(part1_2,mass2),],
                          cpd_split2: [(part2_1,mass1),(part2_2,mass2),],
                        }, 
                    } 
    '''

    for cpd in cpdL:
        #Each split mode
        cpdStr = '+'.join(cpd)
        #Remove redundant split part
        cpd = set(cpd)
        cpdMassL = []
        for part in cpd:
            partL = part.split(';')
            len_partL = len(partL)
            substractLL = [['H'], ['H', 'H2O']]
            for substractL in substractLL:
                cpdMassL.append(\
                    getMass(part, partL[:], len_partL, [], substractL, massRecord))
            #name = part+';-H'
            #if name not in massRecord:
            #    substractH = ['H']
            #    substractH.extend(["H2O"]*(len_partL-1))
            #    mass = massCompute(partL, substractH)
            #    if debug:
            #        print >>sys.stderr, part
            #        print >>sys.stderr, mass
            #cpdMassL.append([name, massRecord[name]])
            #substractH2O = ["H", "H2O"]
            #substractH2O.extend(["H2O"]*(len_partL-1))
            #mass = massCompute(partL, substractH2O)
            #cpdMassL.append([part+';-H-H2O', mass])
        #--------END each part--------------------------------
        massCpdD[cpd_key][cpdStr] = cpdMassL
    #-------END one type of split-----------------------------
#----massSimulate2----------------------------------

def getMass(name, partL, len_partL, addL, substractL, massRecord):
    if addL:
        name += ';+' + '+'.join(addL)
        partL.extend(addL)
    if substractL:
        name += ';-' + '-'.join(substractL)
    #--Minus H2O for each reaction
    substractL.extend(["H2O"]*(len_partL-1))
    mass = massRecord.get(name, '0')
    if mass == '0':
        mass = massCompute(partL, substractL)
        massRecord[name] = mass
    return [name, mass]
    
#-----------------------getMass---------------


def fprint(content):
    """ 
    This is a Google style docs.

    Args:
        param1(str): this is the first param
        param2(int, optional): this is a second param
            
    Returns:
        bool: This is a description of what is returned
            
    Raises:
        KeyError: raises an exception))
    """
    print json_dumps(content,indent=1)
#--------------fprint-------------------

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
        metavar="FILEIN", help="A file to specify the compound \
combinations.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    inputD = readInput(file)
    if debug:
        print >>sys.stderr, "\ninputD\n"
        print >>sys.stderr, inputD

    '''
    compoundL = [[mother_nucleus,R1],
                 [mn,R1,R1],[mn,R1,R1,R1], 
                 [mn,R2,R3], [mn,R2,R3,R2,R3], 
                 [mn,R2,R3,R2,R3,R2,R3]
                ]
    '''
    compoundL = combine(inputD)

    if debug:
        print >>sys.stderr, "\ncompoundL\n"
        print >>sys.stderr, compoundL
    
    masscpdD = {}
    '''
    # Every cpd can be split in multiple ways and would be treated
    # as a list.

    massCpdD(dict): {cpd: {
                          whole: [(mode1, mass1), (mode2, mass2), ]
                          cpd_split1: [(part1_1,mass1),(part1_2,mass2),],
                          cpd_split2: [(part2_1,mass1),(part2_2,mass2),],
                        }, 
                    } 
    '''
    for compound in compoundL:
        massSimulate2(compound, masscpdD)
    
    print '''#First line represents the molecular formula.
#The lines begin with one tab represent different split mode. Each split part is separated by '+'.
#The lines begin with two tabs represent each split part and its mass.
    '''
    for cpd, subD in masscpdD.items():
        print cpd
        #print whole
        print "\twhole"
        for mass in subD['whole']:
            print "\t\t%s\t%.4f" % (mass[0], mass[1])
        #--------------
        split_patternL = subD.keys()
        split_patternL.remove('whole')
        split_patternL.sort()
        for split_pattern in split_patternL:
            print "\t%s" % split_pattern
            massL = subD[split_pattern]
            for mass in massL:
                print "\t\t%s\t%.4f" % (mass[0], mass[1])

    #'''
    #massCpdD(dict): {cpd: [
    #                    [part1, mass1], 
    #                    [part2, mass2]
    #                 ]}                    
    #'''
    #for compound in compoundL:
    #    massSimulate(compound, masscpdD)

    #for cpd, massL in masscpdD.items():
    #    print cpd
    #    for mass in massL:
    #        print "\t%s\t%.4f" % (mass[0], mass[1])
    ##-------------------------------------

    #-------------END reading file----------
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


