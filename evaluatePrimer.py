#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2017, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to evaluate primersearch result.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
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
    parser.add_option("-i", "--primersearch-file", dest="filein",
        metavar="FILEIN", help="Primer search result")
    parser.add_option("-j", "--primersearch-file2", dest="filein2",
        metavar="FILEIN", help="Primer search result2")
    parser.add_option("-k", "--input-file", dest="input_fasta",
        metavar="FILEIN", help="Primer search input fasta")
    parser.add_option("-m", "--input-file2", dest="input_fasta2",
        metavar="FILEIN", help="Primer search input fasta")
    parser.add_option("-p", "--primer-seq", dest="primer",
        metavar="FILEIN", help="Primer file given to primer search")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readInPrimerSearch(file, primerDict, species=''):
    '''
    Primer name SRR037890___comp24_c0_seq1@1
    Amplimer 1
        Sequence: SRR037890___comp24_c0_seq1  
        
        TATTTTCCTATGTTGCTACC hits forward strand at 433 with 0 mismatches
        ACTATTAGCTGTAAAGCAAA hits reverse strand at [23] with 0 mismatches
        Amplimer length: 200 bp


    primerDict = {
        primerName : {

            # Record amplicon info
            # amplicon_start, amplicon_end: 0-started, second excluded
            'species': [
                [targetseq1, AmpliconSize, (amplicon_start, amplicon_end), (forward_mismatch, reverse_mismatch)]
                [targetseq2, AmpliconSize, (amplicon_start, amplicon_end), (forward_mismatch, reverse_mismatch)]
            ],

            # Summary amplicon info
            # Compute number of amplicons for each size of each species
            'sta':{
                species1: {
                    200: 1,
                    300, 1
                },
                species2: {
                    200: 2,
                    300, 1
                }
            }
        }
    }
    '''
    
    seqD = {}
    for line in open(file):
        if line.find("Primer name") == 0:
            primer_name = line[12:-1]
            if primer_name not in primerDict:
                primerDict[primer_name] = {}
                primerDict[primer_name]['sta'] = {}
        elif line.find("Sequence:") != -1:
            targetSeq = line.strip().split()[1]
            seqD[targetSeq] = 1
            if not species:
                species = targetSeq.split('___')[0]
            seq_name = targetSeq
            if species not in primerDict[primer_name]:
                primerDict[primer_name][species] =[]
                primerDict[primer_name]['sta'][species] ={}
        elif line.find("hits forward strand at")  != -1:
            lineL = line.strip().split()
            forward_mismatch = int(lineL[-2])
            amplicon_start = int(lineL[-4])-1
        elif line.find("hits reverse strand at") != -1:
            lineL = line.strip().split()
            reverse_mismatch = int(lineL[-2])
            amplicon_end = (-1) * int(lineL[-4][1:-1])+1
        elif line.find("Amplimer length") != -1:
            ampliconSize = int(line.strip().split()[2])
            primerDict[primer_name][species].append([seq_name, ampliconSize, (amplicon_start, amplicon_end), (forward_mismatch, reverse_mismatch)])
            primerDict[primer_name]['sta'][species][ampliconSize] = primerDict[primer_name]['sta'][species].get(ampliconSize,0)+1
        #----------------------------------------------------------          
        #print >>sys.stderr,primerDict
    #------------------END for---------------------
    return seqD
#-------------------------------------------------------

def tabulizePrimerSearch(primerDict, primerSeqD, query_fastaD, target_fastaD, targetSpeNum):
    '''
    '''
    primeridL = primerDict.keys()
    primeridL.sort(key = lambda x: len(primerDict[x]['sta']), reverse=True)
    headerL = ["Primer ID", "Match species", "Match No", "differentAmpliconCnt", "Amplicon size", "Forward primer", "Reverse Primer", "Amplicon sequence"]
    for primer_id in primeridL:
        speD = primerDict[primer_id]
        staD = speD['sta']
        #clone_spe = len(staD)
        #target = speD['Query_spe']
        matchSpeL = speD.keys()
        matchSpeL.remove("sta")
        matchSpeL.remove('Query_spe')
        clone_spe = len(matchSpeL)
        matchSpeL.sort()
        ampliconSizeL = [subspe+':'+','.join([str(i) for i in staD[subspe].keys()]) for subspe in matchSpeL]
        queryAmpliconSize = 'Query_spe'+':'+ '.'.join([str(i) for i in staD['Query_spe'].keys()])
        ampliconSizeL.insert(0, queryAmpliconSize)
        
        differentAmpliconCnt = len(set(ampliconSizeL))
        
        ampliconsize =';'.join(ampliconSizeL)
        #print >>sys.stderr, primerSeqD[primer_id]
        forward_primer, reverse_primer = primerSeqD[primer_id].split()
        
        #----------------------------------------------------------------------
        querySeqL = speD['Query_spe']
        queryampliconL = []
        for eachQueryMatch in querySeqL:
            id = eachQueryMatch[0]
            start,end = eachQueryMatch[2]
            queryampliconL.append(query_fastaD[id][start:end])
        queryAmplicon = 'Query_spe'+':'+','.join(queryampliconL)
        #-----------------------------------------------------------
        ampliconL = [queryAmplicon]
        for otherspe in matchSpeL:
            targetSeqL = speD[otherspe]
            tmpL = []
            for eachQueryMatch in targetSeqL:
                id = eachQueryMatch[0]
                start,end = eachQueryMatch[2]
                tmpL.append(target_fastaD[id][start:end])
            targetAmplicon = otherspe+':'+','.join(tmpL)
            ampliconL.append(targetAmplicon)
            
        #----------------------------------------------------------------------

        outputL = [primer_id, '; '.join(matchSpeL), str(clone_spe)+'/'+str(targetSpeNum), str(differentAmpliconCnt), 
                ampliconsize, forward_primer, reverse_primer, ';'.join(ampliconL)]
        print '\t'.join(outputL)

    #--------------------------------------------    

#--------------------------------------------------------------------------

def readFasta(file, matchD={},gettargetSpeNum=0):
    seqD = {}
    speS = set()
    for line in open(file):
        if line[0] == '>':
            save = 0
            key = line.strip()[1:]
            if gettargetSpeNum:
                speS.add(key.split('___')[0])
            if matchD.get(key, ''):
                save = 1
                seqD[key] = []
        elif save:
            seqD[key].append(line.strip())
    #----------------------------------------
    for key,seqL in seqD.items():
        seqD[key] = ''.join(seqL)
    if gettargetSpeNum:
        return seqD, len(speS)
    return seqD
#------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    file2 = options.filein2

    input_fasta = options.input_fasta
    input_fasta2 = options.input_fasta2

    primer = options.primer

    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    primerDict = {}
    querySeqD = readInPrimerSearch(file, primerDict, species="Query_spe")
    targetSeqD = readInPrimerSearch(file2, primerDict)

    primerSeqD = dict([line.strip().split('\t',1) for line in open(primer)])
    query_fastaD = readFasta(input_fasta, targetSeqD)
    target_fastaD, targetSpeNum = readFasta(input_fasta2, targetSeqD, gettargetSpeNum=1)

    tabulizePrimerSearch(primerDict, primerSeqD, query_fastaD, target_fastaD, targetSpeNum)

    #print primerDict
    #print json_dumps(primerDict,indent=4)
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
    if verbose:
        print >>sys.stderr,            "--Successful %s" % strftime(timeformat, localtime())

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " %         (' '.join(sys.argv), startTime, endTime)
    fh.close()
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------
#

