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
    This program using <water> to perform local alignment for two
    files.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
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
    parser.add_option("-q", "--query-file", dest="query",
        metavar="QUERY", help="A FASTA file with \
any number of sequences. Each sequence will be compared with \
all sequences in <target-file>. Only strings before the \
first blank will be used as sequence names.")
    parser.add_option("-t", "--target-file", dest="target",
        metavar="TARGET", help="A FASTA file with \
any number of sequences. Only strings before the \
first blank will be used as sequence names.")
    parser.add_option("-o", "--output", dest="out",
        help="The output prefix. One directory named \
by the prefix will be created. and one file <prefix>.table \
will be created.")
    parser.add_option("-a", "--to-all", dest="to_all",
        default=1, help="Output results in one-to-all format or \
one-to-one format. Default 1 indicating <one-to-all>. \
Accept 0 or any ogical FALSE to open <one-to-one> format.")
    parser.add_option("-p", "--parse-only", dest="parse_only",
        default=0, help="Default <0>. Accept 1 to perform only parsing.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.query != None, "A filename needed for -q"
    return (options, args)
#--------------------------------------------------------------------
def readFasta(fasta):
    aDict = {}
    key = ''
    for line in open(fasta):
        if line[0] == '>':
            if key:
                seq = ''.join(seq)
                lenseq = len(seq)
                aDict[key] = [seq, lenseq]
            key = line.split()[0][1:]
            seq = []
        else:
            seq.append(line.strip())
    return aDict
#----END readFasta------------------------------------------

def parseWater(key, lenseq, water, targetD, outTableF):
    fh = open(water)
    while 1:
        line = fh.readline()
        if not line:
            return
        while not line.startswith("#======================================="):
            line = fh.readline()
            if not line:
                return
        line = fh.readline() # blank line
        line = fh.readline() # Aligned_sequences: 2
        line = fh.readline() # # 1: mmu-let-7g-5p
        query_n = line.split()[-1]
        line = fh.readline() # # 2: rev.mmu-let-7g-5p
        target_n = line.split()[-1]
        n = 4
        while n:
            line = fh.readline()
            n -= 1
        #--------Skip 4 lines-----------
        matchLen = fh.readline().split()[2]
        identL = fh.readline().split()
        ident = identL[2].split('/')[0]
        identp = identL[3].strip('()%')

        line = fh.readline() # # Similarity: 
        line = fh.readline() # # Gaps: 
        score = fh.readline().split()[2]
        #print >>outTableF, \
        #    "%s" % '\t'.join([query_n, target_n, str(lenseq), 
        #        matchLen, ident, identp, score])
        # Skip another #==========
        while not line.startswith("#======================================="):
            line = fh.readline()
            if not line:
                return
        #------------------------------------
        line = fh.readline()
        line = fh.readline()
        start = int(line.split()[1])
        #print >>sys.stderr, start, line
        if start <= 3:
            score = str((-1) * float(score))
        print >>outTableF, \
            "%s" % '\t'.join([query_n, target_n, str(lenseq), 
                matchLen, ident, identp, score])

    #------------------------------------
    fh.close()

#----END parseWater----------------------------------------

def executeWater(key, seq, target, targetD, out, outTableF, parse_only):
    output = key+'.water'
    lenseq = len(seq)

    if parse_only:
        parseWater(key, lenseq, out+'/'+output, targetD, outTableF)
        return 
    #-------------------------------------------
    fh = open(key+'.fa', 'w')
    print >>fh, ">%s\n%s" % (key, seq)
    fh.close()
    #----------------------------------------
    if os.path.exists(out+'/'+output):
        return
    cmd = ["water -asequence", key+'.fa', "-bsequence", target, 
       "-gapopen 10 -gapextend 0.5 ", "-outfile", output]
    os.system(' '.join(cmd))
    os.system('/bin/rm -f '+key+'.fa')
    parseWater(key, lenseq, output, targetD, outTableF)
    os.system("/bin/mv " + output + " " + out)
#----------------------------------------------

def executeWaterSep(key, seq, target, targetD, out, outTableF, parse_only):
    if not parse_only:
        fh = open(key+'.fa', 'w')
        print >>fh, ">%s\n%s" % (key, seq)
        fh.close()
    #---------------------------------
    lenseq = len(seq)
    for targetN, targetSL in targetD.items():
        output = key+'-'+targetN+'.water'
        if parse_only:
            parseWater(key, lenseq, out+'/'+output, targetD, outTableF)
        else:
            if os.path.exists(out+'/'+output):
                continue
            tmpfh = open(targetN + '.fa', 'w')
            print >>tmpfh, ">%s\n%s" % (targetN, targetSL[0])
            tmpfh.close()
            cmd = ["water -asequence", key+'.fa', "-bsequence",
                    targetN+'.fa', 
                    "-gapopen 10 -gapextend 0.5 ", "-outfile", output]
            os.system(' '.join(cmd))
            os.system('/bin/rm -f '+targetN+'.fa')
            parseWater(key, lenseq, output, targetD, outTableF)
            os.system("/bin/mv " + output + " " + out)
    os.system('/bin/rm -f '+key+'.fa')
    #----------------------------------------
    # cmd = ["water -asequence", key+'.fa', "-bsequence", target, 
    #    "-gapopen 10 -gapextend 0.5 ", "-outfile", output]
    # os.system(' '.join(cmd))
    # os.system('/bin/rm -f '+key+'.fa')
    # parseWater(key, lenseq, output, targetD, outTableF)
    # os.system("/bin/mv " + output + " " + out)
#------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    query = options.query
    target = options.target
    targetD = readFasta(target)
    out = options.out
    to_all = int(options.to_all)
    os.system("mkdir -p " + out)
    outTable = out+'.table'
    outTableF = open(outTable, 'w')
    parse_only = int(options.parse_only)
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    key = ''
    print >>outTableF, "\t".join(["query", "target", "query_len",
        "alignment_len", "identical_bases", "identity", "score"])
    for line in open(query):
        if line[0] == '>':
            if key:
                if to_all:
                    executeWater(key, ''.join(seq), target, targetD, out,
                        outTableF, parse_only)
                else:
                    executeWaterSep(key, ''.join(seq), target, targetD, out,
                        outTableF, parse_only)
            #-------------------------------------------------
            key = line.split()[0][1:]
            seq = []
        else:
            seq.append(line.strip())
    #----------The last one-------------------- 
    if key:
        if to_all:
            executeWater(key, ''.join(seq), target, targetD, out,
                outTableF, parse_only)
        else:
            executeWaterSep(key, ''.join(seq), target, targetD, out,
                outTableF, parse_only)
    #----------------------------------------------
    outTableF.close()
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


