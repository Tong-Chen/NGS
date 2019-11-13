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
    This is designed to count abnormal chromosomes.

Input:
Label	Result
c2015001	47,XY,+22
c2015002	46,XY
c2015003	47,XX,+16
c2015004	46,XY
c2015005	47,XY,+20
c2015006	69,XXX
c2015007	46,XX
c2015008	47,XY,+4
c2015009	47,XY+16
c2015010	46,XX



'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import re
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

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
        metavar="FILEIN", help="See below in scripts")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def saveD(aDict, key, value):
    if key not in aDict:
        aDict[key] = [value]
    else:
        aDict[key].append(value)
#---------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    aDict = {"Multiploid":{}, "Heteroploid":{}, 
             "Normal":{"Male":[], "Female":[]}, 
             "SingleDiploid":{}, "Microdupdel":{}, 
             "Haploid":{}, "CMultiploid":{}}
    header = 1
    for line in fh:
        if header:
            header -= 1
            continue
        line = line.strip()
        if not len(line): continue # blank
        try:
            if line.find('\t') == -1: continue #only label
        except UnicodeDecodeError, e:
            print >>sys.stderr, e.message
            print >>sys.stderr, line
            sys.exit(1)
        try:
            label, result = line.strip().split('\t')
        except ValueError:
            print >>sys.stderr, line
        if result.find(',') == -1:
            if result == "CMultiploid":
                print "{}\t{}\t{}".format(label, 1, "CMultiploid")
                saveD(aDict["Multiploid"], "CMultiploid", label)
            # Support Tetraploid, but maintainig original construction
            # 2017-09-14, maybe changing later
            elif result == "Tetraploid":
                print "{}\t{}\t{}".format(label, 1, "Tetraploid")
                saveD(aDict["Multiploid"], "Tetraploid", label)
            elif result == "Multiploid":
                # Commented out 20170914
                #print "{}\t{}\t{}".format(label, 1, "Tetraploid")
                #saveD(aDict[result], "Tetraploid", label)
                print "{}\t{}\t{}".format(label, 1, result)
                saveD(aDict[result], result, label)
            elif result == "SingleDiploid":
                print "{}\t{}\t{}".format(label, 4, "Full_CNLOH")
                saveD(aDict[result], "Full_CNLOH", label)
            elif result == "Haploid":
                print "{}\t{}\t{}".format(label, 0, result)
                saveD(aDict[result], result, label)
            elif result.find("CNLOH") != -1:
                print "{}\t{}\t{}".format(label, 4, "Full_CNLOH")
                saveD(aDict["SingleDiploid"], "Full_CNLOH", label)
            else:
                print >>sys.stderr, "##SKIP", line 
            continue
        resultL = result.split(',')
        len_resultL = len(resultL)
        try:
            chrCnt = int(resultL[0])
        except ValueError, e:
            print >>sys.stderr, e.message
            print >>sys.stderr, line
            sys.exit(1)
        sexChr = resultL[1].upper()
        len_sexChr = len(sexChr)
        if chrCnt == 46 and len_sexChr == 2:
            if len_resultL == 2: #c2015178    46,XY
                assert sexChr in ['XX', 'XY'], line
                key = "Female" if sexChr == "XX" else "Male"
                aDict["Normal"][key].append(label)
                print "{}\t{}\t{}".format(label, 5, key)
            else: 
                # Add check of the third element
                # c2015459 46,XY,+14,-21
                # C2016150    49,XX,+5,+7,+11
                #if len_resultL % 2 == 0 and resultL[2][0] in ['+', '-']:  
                # Cancel even number check for C2016150
                if resultL[2][0] in ['+', '-']:  
                    typeL = []
                    for third in resultL[2:]:
                        #commented out in 2017-09-14
                        #assert third[0] in ['+', '-'], line
                        typeL.append(third) 
                    type = ';'.join(typeL)
                    saveD(aDict["Heteroploid"], type, label)
                    print "{}\t{}\t{}".format(label, 2, type)
                else:
                    #commented out in 2017-09-14
                    #assert len_resultL == 3, line

                    #Add loop other mutations 
                    #C2016029 46,XN,dup(3)(q25.31q29)42.60Mb,del(13)(q32.3q34)(15.5Mb)
                    #for third in resultL[2:]:
                    #third = third.upper()
                    #commented out in 2017-09-14
                    third = ','.join(resultL[2:]).upper()
                    if third.find('LOH') != -1: #c2015185    46,XY,CNLOH(14)
                        if third.find('LOH(') != -1:
                            saveD(aDict["SingleDiploid"], "Partial_CNLOH", label)
                            print "{}\t{}\t{}".format(label, 4, "Partial_CNLOH")
                        else:
                            saveD(aDict["SingleDiploid"], "Full_CNLOH", label)
                            print "{}\t{}\t{}".format(label, 4, "Full_CNLOH")
                        #------------SingleDiploid-Full-Partial-----
                    #---Microdupdel----c2015448  46,XX,dup(8)(q22.1q24.3))
                    else: 
                        dup1 = third.find('DUP')
                        del1 = third.find('DEL')
                        if dup1 != -1 and del1 != -1:
                            saveD(aDict['Microdupdel'], 'dupdel', label)
                            print "{}\t{}\t{}".format(label, 3, 'dupdel')
                        elif dup1 != -1:
                            saveD(aDict['Microdupdel'], 'dup', label)
                            print "{}\t{}\t{}".format(label, 3, 'dup')
                        elif del1 != -1:
                            saveD(aDict['Microdupdel'], 'del', label)
                            print "{}\t{}\t{}".format(label, 3, 'del')
                        else:
                            assert 1==0, "Unknown format:" + line
                    #------------------------------------------
                #--------------------------
            #--END chrCnt 46 len_resultL != 2-----
        #------------------------------------------------------------
        elif chrCnt == 46 and len_sexChr != 2: #c1 46,X,+21  or C2 46,XXX,-21
            assert len_resultL == 3, line
            third = resultL[2]
            type = sexChr + third
            saveD(aDict["Heteroploid"], type, label)
            print "{}\t{}\t{}".format(label, 2, type)
        elif chrCnt < 46: #45,XY,-21
            sex_monomer = auto_monomer = 0
            if len(sexChr) == 1:
                sex_monomer = 1
                #saveD(aDict["Heteroploid"], "sex_monomer", label)
            if len_resultL > 2:
                for i in resultL[2:]:
                    i = i.upper()
                    if i[0] == '-':
                        auto_monomer = 1
            if sex_monomer and auto_monomer:
                saveD(aDict["Heteroploid"], "sex_auto_monomer", label)
                print "{}\t{}\t{}".format(label, 2, "sex_auto_monomer")
            elif sex_monomer:
                saveD(aDict["Heteroploid"], "sex_monomer", label)
                print "{}\t{}\t{}".format(label, 2, "sex_monomer")
            elif auto_monomer:
                saveD(aDict["Heteroploid"], "auto_monomer", label)
                print "{}\t{}\t{}".format(label, 2, "auto_monomer")
            else:
                assert 1==0, "Unknown format: " + line
            #----------------------------------------------------
        elif chrCnt < 65:  #48,XXY,+21 #48,XY,+12,+18
            trisom_cnt = chrCnt - 46
            trisom_typL = []
            if trisom_cnt not in aDict["Heteroploid"]:
                aDict["Heteroploid"][trisom_cnt] = {}
            if len_sexChr > 2:
                trisom_typL.append(sexChr)
            for trism_auto in resultL[2:]:
                trism_auto = trism_auto.strip().upper()
                if trism_auto[0] == '+':
                    trisom_typL.append(trism_auto[1:])
                else:
                    assert 'DUP' in trism_auto or 'DEL' in trism_auto, line
            trisom_typL.sort()
            trisom_typ = ';'.join(trisom_typL)
            saveD(aDict["Heteroploid"][trisom_cnt], trisom_typ, label)    
            print "{}\t{}\t{}:{}".format(label, 2, trisom_cnt, trisom_typ)
        elif chrCnt < 69 or chrCnt > 69:
            n3 = chrCnt-69
            if n3 > 0:
                type = '3N+'+str(n3)
            else:
                type = '3N'+str(n3)
            if type not in aDict["Multiploid"]:
                aDict["Multiploid"][type] = {}
            mp_typL = [sexChr]
            for mp_type in resultL[2:]:
                mp_type = mp_type.strip().upper()
                # skip line like C2016302  70,XXY,+10,DUP(X)(p22.32p22.13)(13.5Mb))))
                if mp_type.find('DUP') != -1 or mp_type.find('DEL') != -1 \
                        or mp_type.find('CNL') != -1:
                            continue
                if n3 > 0:
                    assert mp_type[0] == '+', line
                else:
                    assert mp_type[0] == '-', line
                mp_typL.append(mp_type)
            mp_typ = ';'.join(mp_typL)
            saveD(aDict["Multiploid"][type], mp_typ, label)
            print "{}\t{}\t{}".format(label, 1, type+':'+mp_typ)
        elif chrCnt == 69:
            # Commented out at 20170914 for
            # C2016637  69, XXY, CNLOH(X)(p22.2q22.1)(61Mb)
            #assert len_resultL == 2, line
            assert len_resultL >= 2, line
            saveD(aDict["Multiploid"], str(chrCnt)+','+sexChr, label)
            print "{}\t{}\t{}:{}".format(label, 1, chrCnt, sexChr)
        else:
            assert 1==0, "Unknown format" + line

    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    orderL = ["Multiploid", "Heteroploid", 
              "Microdupdel", "Haploid", 
              "SingleDiploid", "Normal"]
    countD = {}
    #print aDict
    for type, subD in aDict.items():
        if type not in countD:
            countD[type] = {}
            countD[type]['total'] = 0
        for subType, subVL in subD.items():
            countD[type][subType] = {}
            countD[type][subType]['total'] = 0
            if isinstance(subVL, list):
                count = len(subVL)
                countD[type][subType]['total'] += count
                countD[type]['total'] += count
            elif isinstance(subVL, dict):
                for subsubType, subsubVL in subVL.items():
                    assert isinstance(subsubVL, list), subsubVL
                    count = len(subsubVL)
                    countD[type][subType][subsubType] = count
                    countD[type][subType]['total'] += count
                    countD[type]['total'] += count

    for type in orderL:
        print "{}\t{}".format(type, countD[type]['total'])
        for subType in aDict[type].keys():
            print "\t{}\t{}".format(subType, countD[type][subType]['total'])
            subVD = aDict[type][subType]
            if isinstance(subVD, dict):
                for subsubType in subVD.keys():
                    print "\t\t{}\t{}".format(subsubType, countD[type][subType][subsubType])

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

'''
Label	Result
c2015001	47,XY,+22
c2015002	46,XY
c2015003	47,XX,+16
c2015004	46,XY
c2015005	47,XY,+20
c2015006	69,XXX
c2015007	46,XX
c2015008	47,XY,+4
c2015009	47,XY+16
c2015010	46,XX
c2015011	46,XX
c2015012	47,XY,+22
c2015013	47,XY,+16
c2015014	46,XX
c2015015	46,XX
c2015016	47,XY,+4
c2015017	46,XY
c2015018	47,XX,+22
c2015019	47,XX,+16
c2015020	46,XX
c2015021	46,XX
c2015022	45,X
c2015023	46,XY
c2015024	47,XX,+16
c2015025	46,XX
c2015026	46,XX
c2015027	47,XY,+7
c2015028	46,XY
c2015029	47,XY,+16
c2015030	47,XX,+22,dup(9)(p23)2.9Mb
c2015031	47,XY,+10,del(6)(q25.3q27)(22Mb)dup(7)(q31.32q36.3)159Mb
c2015032	47,XX,+15
c2015033	47,XY,+14
c2015034	47,XX,+22
c2015035	47,XXX
c2015036	47,XX,+8
c2015037	47,XX,+22
c2015038	46,XY
c2015039	47,XY,+14
c2015040	46,XY
c2015041	69,XXY
c2015042	46,XY,dup(14)(q12)(2.4Mb)
c2015043	47,XX,+22
c2015044	46,XX
c2015045	46,XY
c2015046	47,XY+18
c2015047	46,XY
c2015048	46,XY
c2015049	47,XX,+15
c2015050	45,X
c2015051	46,XY
c2015052	46,XX
c2015053	Multiploid
c2015054	47,XX,+9
c2015055	46,XX
c2015056	46,XY
c2015057	46,XY
c2015058	47,XX,+16
c2015059	47,XY,+16
c2015060	47,XY,+22
c2015061	47,XY,+22
c2015062	47,XX,+16
c2015063	45,X
c2015064	Multiploid
c2015065	46,XX
c2015066	47,XY,+15
c2015067	47,XY,+16
c2015068	46,XX
c2015069	47,XX,+16
c2015070	48,XY,+15,+18
c2015071	47,XX,+22
c2015072	46,XY
c2015073	47,XY,+20
c2015074	69,XXY
c2015075	47,XY,+15
c2015076	47,XY,+16
c2015077	45,X
c2015078	46,XX,del(8)(p23.3p21.2)26Mbdup(8)(p21.2q24.3)119Mbdup(20)(p13p12.3)5.8Mbdup(20)(p12.1)4.2Mb
c2015079	47,XY,+5
c2015080	46,XX
c2015081	47,XY,+15
c2015082	47,XY,+16
c2015083	46,XY
c2015084	47,XY,+10
c2015085	Multiploid
c2015086	46,XX
c2015087	46,XX
c2015088	47,XY,+15
c2015089	47,XX,+14
c2015090	47,XX,+21
c2015091	46,XX
c2015092	46,XY
c2015093	46,XY
c2015094	46,XY
c2015095	47,XY,+9
c2015096	46,XX
c2015097	46,XX
c2015098	46,XX
c2015099	47,XX,+16
c2015100	47,XY,+16
c2015101	46,XY
c2015102	46,XX
c2015103	46,XX
c2015104	46,XX
c2015105	48,XX,+7,+16
c2015106	47,XX,+21
c2015107	47,XX,+21
c2015108	47,XY,+16
c2015109	46,XX
c2015110	47,XX,+15
c2015111	47,XX,+22
c2015112	69,XXX
c2015113	46,XY
c2015114	47,XY,+2
c2015115	46,XY
c2015116	47,XX,+22
c2015117	47,XX,+16
c2015118	48,XY,+18,+22
c2015119	46,XY
c2015120	46,XY
c2015121	45,X
c2015122	48,XY,+16,+22
c2015123	47,XX,+22
c2015124	47,XY,+4
c2015125	46,XX
c2015126	46,XY
c2015127	46,XX
c2015128	47,XX,+15
c2015129	47,XX,+5
c2015130	47,XX,+5
c2015131	46,XX
c2015132	46,XX,CNLOH(2)(q14.1q24.3)43MbCNLOH(X)(q23q28)31Mb
c2015133	46,XY
c2015134	47,XX,+16
c2015135	46,XX
c2015136	45,X
c2015137	47,XY,+16
c2015138	47,XX,+21
c2015139	46,XY,del(1)(p32.3)2.5Mb
c2015140	47,XX,+16
c2015141	47,XY,+16
c2015142	46,XX
c2015143	Multiploid
c2015144	47,XY,+13
c2015145	46,XY
c2015146	46,XY
c2015147	47,XX,+16
c2015148	47,XX,+16
c2015149	47,XY,+8
c2015150	46,XX
c2015151	47,XX,+4
c2015152	47,XX,+2
c2015153	46,XX
c2015154	46,XX
c2015155	46,XY
c2015156	48,XY,+15,+21
c2015157	70,XXX,+11
c2015158	46,XY,del(18)(p11.32p11.21)13Mb
c2015159	48,XY,+21x2
c2015160	46,XY
c2015161	47,XY,+16
c2015162	48,XY,+15,+16
c2015163	47,XY,+9
c2015164	Multiploid
c2015165	45,X
c2015166	46,XY
c2015167	45,X
c2015168	46,XX
c2015169	47,XY,+16
c2015170	Multiploid
c2015171	46,XY
c2015172	46,XX
c2015173	45,X
c2015174	47,XX,+14
c2015175	46,XX
c2015176	Multiploid
c2015177	48,XX,+2,+15
c2015178	46,XY
c2015179	46,XX
c2015180	46,XX
c2015181	47,XY,+16
c2015182	46,XX
c2015183	47,XX,+15
c2015184	68,XXY,-9
c2015185	46,XY,CNLOH(14)
c2015186	46,XY
c2015187	46,XY
c2015188	46,XX
c2015189	
c2015190	46,XY
c2015191	Multiploid
c2015192	47,XX,+18
c2015193	45,X
c2015194	46,XX
c2015195	47,XX,+22
c2015196	47,XX,+21
c2015197	47,XX,+22
c2015198	46,X,+22
c2015199	46,XX
c2015200	46,XX
c2015201	46,XX
c2015202	47,XX,+13
c2015203	47,XY,+4
c2015204	46,XX
c2015205	46,XX
c2015206	46,XX
c2015207	46,XX
c2015208	46,XX
c2015209	47,XY,+22
c2015210	46,XX
c2015211	46,XY
c2015212	46,XX
c2015213	46,XY
c2015214	47,XX,+13
c2015215	47,XY,+22
c2015216	47,XX,+22
c2015217	46,XX,del(15)(q26.1q26.3)(8.1Mb)
c2015218	47,XX,+16
c2015219	47,XY,+16
c2015220	46,XX
c2015221	45,X
c2015222	46,XX
c2015223	45,X
c2015224	47,XX,+22
c2015225	47,XX,+16
c2015226	47,XY,+21
c2015227	46,XY
c2015228	47,XY,+17
c2015229	47,XX,+16
c2015230	47,XY,+16
c2015231	46,XX
c2015232	47,XX,+21
c2015233	46,XY
c2015234	46,XY
c2015235	Multiploid
c2015236	47,XY,+16
c2015237	47,XY,+22
c2015238	45,X
c2015239	46,XX
c2015240	46,X,+21
c2015241	46,XY
c2015242	Multiploid
c2015243	47,XY,+16
c2015244	69,XXX
c2015245	46,XX
c2015246	47,XX,+3
c2015247	46,XX
c2015248	47,XY,+15
c2015249	47,XY,+14
c2015250	46,XX
c2015251	47,XX,+7
c2015252	47,XY,+13
c2015253	47,XY,+22
c2015254	46,XY
c2015255	46,XY
c2015256	46,XX
c2015257	47,XX,+21
c2015258	47,XX,+16
c2015259	46,XX
c2015260	46,XY,+9
c2015261	69,XXX
c2015262	45,X
c2015263	46,XY
c2015264	47,XX,+13
c2015265	47,XX,+15
c2015266	47,XY,+4
c2015267	46,XY
c2015268	47,XX,+16
c2015269	69,XXY
c2015270	46,XX
c2015271	69,XXY
c2015272	47,XY,+16
c2015273	46,XX
c2015274	47,XX,+21
c2015275	47,X,del(X)(q26.2q28)23.7Mb,+16
c2015276	45,X
c2015277	46,XX
c2015278	46,XX
c2015279	46,XX
c2015280	46,XY,dup(8)(q21.2q24.23)58Mb
c2015281	46,XX
c2015282	46,XX
c2015283	47,XX,+17
c2015284	47,XY,+16
c2015285	47,XY,+16
c2015286	47,XY,+16
c2015287	46,X,+13
c2015288	46,XY
c2015289	46,XY
c2015290	47,XY,+21
c2015291	69,XXY
c2015292	Multiploid
c2015293	47,XY,+22
c2015294	47,XX,+15
c2015295	47,XY,+8
c2015296	46,XY,dup(8)(q24.13q24.3)20.1Mbdel(14)(q32.2q32.33)6.3Mb
c2015297	46,XX
c2015298	47,XY,+22
c2015299	46,XX,UPD(5)
c2015300	47,XX,+9
c2015301	46,XX
c2015302	47,XY,+16
c2015303	46,XY
c2015304	46,XY
c2015305	47,XX,+14
c2015306	47,XX,+16
c2015307	46,XY
c2015308	46,XY
c2015309	46,XY
c2015310	48,XXX,+14
c2015311	46,XX
c2015312	46,XY
c2015313	46,XY
c2015314	46,XY
c2015315	69,XXY
c2015316	46,XY
c2015317	46,XX,dup(2)(p25.3p25.2)3.6Mb
c2015318	46,XX
c2015319	46,XY
c2015320	46,XY
c2015321	46,XX
c2015322	46,XY,del(12)(q12q13.11)9Mb
c2015323	48,XY,+7,+11
c2015324	46,XX
c2015325	47,XX,+16
c2015326	46,XX,del(18)(p11.32p11.21)14Mb
c2015327	47,XY,+4
c2015328	47,XX,+22
c2015329	47,XX,+13
c2015330	46,XY
c2015331	46,XX
c2015332	47,XY,UPD(4)58Mb,+7,UPD(9)(q21.11q31.3)40Mb
c2015333	47,XY,+12
c2015334	46,XY
c2015335	46,XX
c2015336	47,XX,+15
c2015337	46,XY
c2015338	46,XX
c2015339	46,XY
c2015340	46,XX
c2015341	47,XX,+16
c2015342	46,XX
c2015343	46,XX
c2015344	46,XY
c2015345	46,XX
c2015346	46,XY
c2015347	47,XX,+22
c2015348	46,XY
c2015349	46,XX
c2015350	47,XY,+7
c2015351	46,XX
c2015352	47,XY,+16
c2015353	47,XY,+16
c2015354	46,XX
c2015355	47,XY,+17
c2015356	47,XY,+20
c2015357	46,XY
c2015358	46,XX
c2015359	47,XY,+22
c2015360	46,XY
c2015361	46,XY
c2015362	46,XX
c2015363	47,XX,+2
c2015364	48,XX,+21×2
c2015365	47,XX,+9
c2015366	47,XY,+22
c2015367	48,XXX,+2,dup(7)(q11.22136.3)90Mb
c2015368	47,XX,+16
c2015369	46,XX
c2015370	46,XX
c2015371	46,XY
c2015372	47,XX,+13
c2015373	46,XX
c2015374	46,XY
c2015375	47,XY,+16
c2015376	47,XY,+22
c2015377	47,XX,+22
c2015378	69,XXY
c2015379	46,XX
c2015380	47,XY,+16
c2015381	46,XY
c2015382	46,XX
c2015383	48,XY,+11,+22
c2015384	47,XX,+15
c2015385	47,XX,+16
c2015386	47,XY,+22
c2015387	47,XY,+3
c2015388	47,XY,+12
c2015389	46,XY
c2015390	46,XY
c2015391	47,XY,+16
c2015392	46,XX
c2015393	46,XX
c2015394	47,XY,+10
c2015395	46,XY
c2015396	45,X
c2015397	46,XY
c2015398	Multiploid
c2015399	XY,del(8)(p23.3p22)2.1Mbdup(8)(p23.2)2.2Mbdup(18)(q12.3q23)40.7Mb
c2015400	47,XY,+21
c2015401	47,XX,+5
c2015402	46,XY
c2015403	47,XX,+22
c2015404	47,XY,+16
c2015405	46,XX,del(2)(q37.1q37.3)9.2Mbdup(6)(p25.3p24.2)11.3Mb
c2015406	47.XX,+18
c2015407	46,XY
c2015408	Multiploid
c2015409	48,XY,+5,+20
c2015410	47,XY,+4
c2015411	46,XX
c2015412	46,XX
c2015413	47,XX,+16
c2015414	46,XX
c2015415	47,XY,+16
c2015416	47,XX,+10
c2015417	47,XY,+16
c2015418	46,XX
c2015419	45,XY,-21
c2015420	47,XX,+16
c2015421	48,XXX,+6
c2015422	46,XX
c2015423	47,XY,+16
c2015424	47,XY,+16
c2015425	47,XY,+16
c2015426	46,XY
c2015427	Multiploid
c2015428	46,XY
c2015429	47,XY,+9
c2015430	47,XY,+21
c2015431	47,XY,+23
c2015432	47,XX,+12
c2015433	47,XY,+16
c2015434	46,XX
c2015435	46,XX
c2015436	47,XX,+22
c2015437	46,XY
c2015438	48,XXY,+21
c2015439	48,XY,+12,+18
c2015440	46,XX
c2015441	46,XX
c2015442	47,XY,+15
c2015443	47,XY,+16
c2015444	46,XX
c2015445	46,XY
c2015446	47,XX,+21
c2015447	46,XX
c2015448	46,XX,dup(8)(q22.1q24.3)42.4Mb,dup(13),(q31.1q33.3)22.3Mbdel(13)(q33.3q34)16.8Mb
c2015449	46,XX
c2015450	45,XO,dup(22)(q11.21)3.3Mb
c2015451	46,XX,del(3)(q24q26.32)32Mb
c2015452	46,XX
c2015453	47,XX,+13
c2015454	46,XX
c2015455	46,XX
c2015456	46,XX
c2015457	47,XX,+13
c2015458	47,XX,+16
c2015459	46,XY,+14,-21
c2015460	46,XY,del(4)(p16.2)3.7Mbdup(4)(p16.3p15.1)28.4Mb
c2015461	47,XX,+16
c2015462	47,XY,+16
c2015463	47,XX,+9
c2015464	48,XY,dup(3)(q26.1q29)30.6Mb,+18,+21
c2015465	47,XX,+20
c2015466	46,XY
c2015467	46,XX
c2015468	46,XY
c2015469	46,XY
c2015470	46,XY
c2015471	47,XY,+8
c2015472	48,XX,+7,+20
c2015473	46,XX
c2015474	47,XY,+16
c2015475	46,XX
c2015476	69,XXX
c2015477	47,XX,+22
c2015478	47,XY,+15
c2015479	45,X
c2015480	46,XY
c2015481	47,XX,+22
c2015482	46,XY
c2015483	46,XX
c2015484	47,XX,+18
c2015485	46,XX
c2015486	46,XX
c2015487	48,XY,+6,+18
c2015488	69,XXY
c2015489	68,XXY,-4
c2015490	46,XY
c2015491	46,XX
c2015492	47,XY,+13
c2015493	47,XX,+22
c2015494	47,XX,+2
c2015495	46,XY
c2015496	46,XX
c2015497	47,XY,+22
c2015498	46,XX,dup(10)(q25.3q26.3)18.1Mb
c2015499	47,XX,+16
c2015500	46,XX
c2015501	47,XX,+22
c2015502	
c2015503	46,XX
c2015504	47,XX,+16
c2015505	47,XY,+22
c2015506	47,XX,+19
c2015507	69,XXX
c2015508	46,XX
c2015509	46,XY
c2015510	47,XX,+16
c2015511	46,XX,del(8)(p23.3p21.3)20Mbdup(21)(q21.3q22.3)22.9Mb
c2015512	46,XY
c2015513	47,XX,+14,del(9)(q34.2q34.3)3.76Mbdup(9)(34.13q34.2)1.8Mb
c2015514	Multiploid
c2015515	47,XY,+15
c2015516	46,XY
c2015517	91,XXX
c2015518	46,XY
c2015519	48,+16×2
c2015520	46,XX
c2015521	45,X
c2015522	46,XX
c2015523	46,XY
c2015524	Multiploid
c2015525	46,XY
c2015526	46,XX
c2015527	47,XY,+12
c2015528	47,XY,+16
c2015529	46,XY
c2015530	46,XY
c2015531	47,XY,+22
c2015532	46,XX
c2015533	47,XY,+22
c2015534	46,XY
c2015535	47,XX,+21
c2015536	46,XY
c2015537	46,XY
c2015538	46,XX
c2015539	46,XX
c2015540	47,XX,+8
c2015541	46,XX
c2015542	47,XX,+20
c2015543	49,XY,+16×3
c2015544	Multiploid
c2015545	47,XX,+22
c2015546	46,XX
c2015547	47,XY,+22
c2015548	69,XXX
c2015549	46,XY
c2015550	46,XX
c2015551	46,XY
c2015552	47,XY,+20
c2015553	47,XX,+21
c2015554	Multiploid
c2015555	46,XX
c2015556	46,XX
c2015557	47,XX,+6
c2015558	46,XY
c2015559	46,XY
c2015560	46,XX,dup(3)(p26.3,p14.3)56.7Mbdel(8)(p23.3p21.27)22Mb
c2015561	46,XY
c2015562	46,XX
c2015563	Multiploid
c2015564	46,XY
c2015565	46,XX
c2015566	47,XY,+16
c2015567	46,XX
c2015568	Multiploid
c2015569	46,XX
c2015570	45,X
c2015571	46,XX
c2015572	47,XX,+16
c2015573	47,XY,+16
c2015574	46,XY
c2015575	46,XX
c2015576	46,XX
c2015577	47,XX,+6
c2015578	46,XX
c2015579	45,XO
c2015580	47,XXX
c2015581	47,XY,+7
c2015582	47,XY,+9
c2015583	46,XY
c2015584	46,XX
c2015585	47,XX,+2
c2015586	46,XX
c2015587	47,XY,+11
c2015588	47,XX,+4
c2015589	46,XY,del(18)(p11.1pten?)
c2015590	47,XX,+16
c2015591	49,XX,+16×3
c2015592	46,XX
c2015593	46,XX,del(21)(q11.2q21.1)8.2Mbdup(10)(pten?P21.1)56Mb
c2015594	46,XX
c2015595	69,XXY
c2015596	46,XX
c2015597	46,XX
c2015598	46,XX
c2015599	47,XX,+3
c2015600	47,XX,+22
c2015601	46,XY
c2015602	46,XY
c2015603	46,XX
c2015604	69,XXX
c2015605	48,XX,+13,+16
c2015606	47,XY,+15
c2015607	47,XY,+16
c2015608	46,XX
c2015609	48,XY,+14,+20
c2015610	46,XY
c2015611	47,XY,+16
c2015612	48,XY,+15,+18
c2015613	47,XY,+4
c2015614	47,XY,+12
c2015615	47,XX,+16
c2015616	46,XX,del(7)(q11.23q21.11)9.3Mb
c2015617	47,XX,+16
c2015618	46,XY
c2015619	46,XX,dup(17)(q12)1.5Mb
c2015620	47,XY,+2
c2015621	47,XX,+22
c2015622	46,XX
c2015623	47,XXX
c2015624	46,XX
c2015625	46,XY
c2015626	45,XO
c2015627	47,XX,+16
c2015628	46,dup(2)(q33.3q37.3)31Mbdel(2)(q37.3)301Mb
c2015629	46,XX,+14,-21
c2015630	45,XO
c2015631	46,XX
c2015632	48,XX,+15,+16
c2015633	46,XY
c2015634	46,XX
c2015635	46,XY
c2015636	47,XX,+16
c2015637	50,XYY,+12,+15,+20
c2015638	47,XX,+4
c2015639	48,XX,+16,+20
c2015640	46,XX
c2015641	Multiploid
c2015642	45,X
c2015643	48,XY,+14,+22
c2015644	69,XXY
c2015645	46,XY
c2015646	47,XX,+2
c2015647	49,+21×3
c2015648	69,XXX
c2015649	46,XY
c2015650	46,XX,dup(2)(q11.1q37.3)47.5Mbdel(10)(q26.2q26.3)7,1Mb
c2015651	46,XY
c2015652	46,XY
c2015653	47,XX,+16
c2015654	46,XX,upd(21)(q21.3q22.11)7.2Mb
c2015655	46,XX
c2015656	45,XX,-21
c2015657	46,XY,del(10)(q26.2q26.3)4.6Mb
c2015658	46,XY
c2015659	47,XY,+16
c2015660	46,XY
c2015661	46,XY
c2015662	46,XX,-18p,+18q
c2015663	46,XX
c2015664	69,XXY
c2015665	46,XY
c2015666	47,XX,+16
c2015667	46,XX
c2015668	47,XY,+16
c2015669	46,XX
c2015670	46,XX,dup(6)(q24q27)30.7Mb
c2015671	47,XY,+7
c2015672	46,XY
c2015673	46,XY,udp(2)
c2015674	47,XX,+7
c2015675	Multiploid
c2015676	47,XY,+22
c2015677	46,XY
c2015678	47,XY,+3
c2015679	46,XY
c2015680	46,XY
c2015681	46,XX
c2015682	47,XY,+16
c2015683	48,XX,+21×2
c2015684	46,XY,母46,XX,t(16,19)(p11.2p12)
c2015685	47,XX,+22
c2015686	46,XX
c2015687	Multiploid
c2015688	46,XY
c2015689	46,XY
c2015690	47,XY,+16
c2015691	47,XY,+22
c2015692	46,XY
c2015693	69,XXX
c2015694	47,XY,+22
c2015695	47,XX,+16
c2015696	46,XY
c2015697	47,XY,+22
c2015698	47,XY,+22
c2015699	46,XX
c2015700	47,XY,+10
c2015701	46,XY,del(6)(q25.2q27)17.3Mbdup(16)(q23.1q24.3)14.1Mb
c2015702	47,XX,+16
c2015703	46,XX
c2015704	48,XY,+5,+6
c2015705	Multiploid
c2015706	48,XXX,+12
c2015707	46,XY
c2015708	46,XX
c2015709	46,XX
c2015710	47,XY,+16
c2015711	46,XY
c2015712	45,X
c2015713	46,XX
c2015714	46,XY
c2015715	47,XX,+16
c2015716	46,XY
c2015717	45,XY,-21
c2015718	46,XY
c2015719	
c2015720	47,XY,+3
c2015721	47,XY,+5
c2015722	47,XY,+13
c2015723	47,XY,+13
c2015724	45,X
c2015725	Multiploid
c2015726	47,XX,+16
c2015727	46,XY
c2015728	47,XY,+22
c2015729	47,XY,+10
c2015730	47,XX,+16
c2015731	46,XX
c2015732	47,XX,+21
c2015733	47,XY,+16
c2015734	46,XY
c2015735	45,X
c2015736	47,XX,+2
c2015737	46,XY
c2015738	46,XX
c2015739	47,XX,+15
c2015740	46,XX
c2015741	47,XY,+14
c2015742	46,XY
c2015743	47,XX,+14
c2015744	47,XX,+15
c2015745	47,XX,+16
c2015746	46,XY,del(6)(p25.3p25.2)2.3Mb
c2015747	46,XX,del(4)(p16.3p15.1)30.7Mb
c2015748	46,XY,del(4)(q34.3q35.2)9.0Mb
c2015749	47,XY,+22
c2015750	47,XY,+16
c2015751	46,XY
c2015752	46,XX
c2015753	47,XX,+3
c2015754	46,XX
c2015755	46,XY
c2015756	46,XX
c2015757	47,XY,dup(5)(p15.33p13.3)30Mb
c2015758	46,XY
c2015759	46,XX
c2015760	45,X
c2015761	45,X
c2015762	69,XXX
c2015763	46,XX
c2015764	46,XX
c2015765	45,X
c2015766	Multiploid
c2015767	46,XX
c2015768	46,XX
c2015769	46,XY
c2015770	47,XY,+22
c2015771	47,XX,+21
c2015772	47,XX,+9
c2015773	46,XY
c2015774	47,XY,+22
c2015775	46,XX
c2015776	46,XX,+12
c2015777	47,XY,+22
c2015778	46,XX
c2015779	46,XX
c2015780	47,XY,+8
c2015781	46,XY
c2015782	45,X/46,XY
c2015783	47,XY,+3
c2015784	47,XY,+16
c2015785	SingleDiploid
c2015786	47,XY,+16
c2015787	46,XY
c2015788	47,XX,+16
c2015789	47,XY,+13
c2015790	49,XY,+7,+18+20
c2015791	46,XY
c2015792	47,XY,+16
c2015793	46,XY
c2015794	46,XY
c2015795	46,XX
c2015796	47,XY,+22
c2015797	46,XY,del(6)(q25.3q27)10.6Mb
c2015798	45,XY,-21
c2015799	47,XY,+8
c2015800	46,XX
c2015801	46,XY
c2015802	46,XX
c2015803	47,XY,+16
c2015804	46,XX
c2015805	47,XX,+22
c2015806	46,XX
c2015807	47,XY,+16
c2015808	46,XX
c2015809	46,XY
c2015810	45,X
c2015811	46,XY
c2015812	47,XY,+21
c2015813	46,XY,del(16)(p13.11p12.3)3.3Mb
c2015814	46,XY
c2015815	46,XY
c2015816	46,XY
c2015817	47,XYY
c2015818	46,XY
c2015819	47,XX,+8
c2015820	47,XX,+16
c2015821	47,XY,+16
c2015822	48,XY,+9,+13
c2015823	Multiploid
c2015824	46,XX
c2015825	47,XX,+10
c2015826	47,XY,+3,dup(x)(q27.1q28)15Mb
c2015827	47,XX,+13
c2015828	47,XX,+22
c2015829	47,XY,+16
c2015830	47,XY,+16
c2015831	46,XY
c2015832	69,XXX
c2015833	46,XY
c2015834	47,XY,+22
c2015835	46,XY
c2015836	47,XY,+21
c2015837	46,XY
c2015838	47,XX,+16
c2015839	46,XY
c2015840	46,XY
c2015841	46,XX
c2015842	46,XY
c2015843	46,XX,dup(15)(q15.1q26.3)59.4Mb
c2015844	46,XX
c2015845	45,X
c2015846	47,XX,+22
c2015847	47,XY,+2,dup(13)(q12.12)13Mb
c2015848	47,XY,+17
c2015849	45,X
c2015850	Multiploid
c2015851	46,XY
c2015852	45,X
c2015853	Multiploid
c2015854	47,XY,+16
c2015855	46,XX
c2015856	47,XY,+13
c2015857	47,XX,+16
c2015858	46,XY
c2015859	46,XY
c2015860	47,XX,+6
c2015861	47,XX,+6
c2015862	69,XXY
c2015863	46,XX
c2015864	46,XY
c2015865	46,XY
c2015866	48,XXY,+16
c2015867	46,XY,del(11)(q32q44)(42AB),dup(11)(q23.3,q25)15.Mb
c2015868	47,XX,+2
c2015869	47,XY,+16
c2015870	47,XX,+16
'''
