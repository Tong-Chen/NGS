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
    This is designed to run psiblast and parse psiblast result.
'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0


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
        metavar="FILEIN", help="Query fasta file")
    parser.add_option("-d", "--target-db", dest="target_db",
        help="Target DB")
    parser.add_option("-m", "--identity", dest="identity",
        default=50, type=float, help="Identity like 50 (default) representing larger than 50% identity")
    parser.add_option("-c", "--coverage", dest="coverage", default=50, 
        type=float, help="Coverage like 50 (default) representing larger than 50% coverage of query seq")
    parser.add_option("-e", "--evalue", dest="evalue",
        default='0.0001', help="A number, 0.0001 (default)")
    parser.add_option("-p", "--parse-only", dest="parse_only",
        default=False, action="store_true", help="Specify to only do parsing")
    parser.add_option("-P", "--python-parse-cmd", dest="parse_cmd",
        default=True, action="store_false", help="Use python instead of shell do parsing.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def parsePsiblastResult(file, identity_thresh, coverage_thresh):
    aDict = {}
    for line in open(file):
        if line[0] == '#':
            continue
        #--------------------------
        line = line.strip()
        if not line:
            break
        lineL = line.split("\t")
        if debug:
            print >>sys.stderr, lineL
        identity_value = float(lineL[5])
        coverage_value = float(lineL[6])
        evalue = lineL[7]
        if identity_value < identity_thresh or coverage_value < coverage_thresh:
            continue
        key = ';'.join(lineL[0:2])
        if key not in aDict:
            aDict[key] = {}
            aDict[key]['identity'] = identity_value
            aDict[key]['evalue'] = evalue
            aDict[key]['coverage'] = coverage_value
            aDict[key]['result'] = lineL[1:]
        #-----------------------------------
        else:
            if identity_value > aDict[key]['identity'] and \
                    coverage_value > aDict[key]['coverage']:
                aDict[key]['identity'] = identity_value
                aDict[key]['evalue'] = evalue
                aDict[key]['coverage'] = coverage_value
                aDict[key]['result'] = lineL[1:]
            #----------------------------
        #---------------------------------
    #---------------------------------
    return aDict
#--------------------------------

def outputPsiResult(aDict, fh=sys.stdout):
    keyL = aDict.keys()
    keyL.sort(key=lambda x: [aDict[x]['evalue'], (-1)*aDict[x]['identity'], (-1)*aDict[x]['coverage']])
    for key in keyL:
        print >>fh, '\t'.join(aDict[key]['result'])
#---------------------------------
def parsePsiblastResultCMD(file, identity_thresh, coverage_thresh, output):
    identity_thresh = str(identity_thresh)
    coverage_thresh = str(coverage_thresh)
    cmd = ["grep '^1'", file, "| sort -k8,8n -k6,6nr -k7,7nr | awk 'BEGIN{OFS=FS=\"\t\"}{if(save[$2]==\"\" && $6>="+identity_thresh+" && $7>="+coverage_thresh+") {print $0; save[$2]=1;}}' | cut -f 1 --complement | sed '1 i\Subject\tQuery length\tSubject length\tAlignment length\tPercentage of identical matches\tQuery Coverage Per Subject\tevalue\tTaxonomyID\tTaxonomy name\tCommon name' >", output]
    if debug:
        print ' '.join(cmd)
    os.system(' '.join(cmd))
#--------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    target_db = options.target_db
    identity  = options.identity
    coverage  = options.coverage
    evalue    = options.evalue
    parse_only = options.parse_only
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    parse_cmd = options.parse_cmd
    if parse_only:
        if parse_cmd:
            parsePsiblastResultCMD(file+'.psiblast.table',identity,coverage, file+'.psiblast.table.parsed.tsv')
            return
        aDict = parsePsiblastResult(file+'.psiblast.table', identity, coverage)
        final = file+'.psiblast.table.parsed.tsv'
        final_fh = open(final, 'w')
        outputPsiResult(aDict, fh=final_fh)
        final_fh.close()
        return
    #-----------------------------------------------------
    fasta_number = int(os.popen("grep -c '>' "+file).read())
    if debug:
        print >>sys.stderr, fasta_number
    if fasta_number > 1:
        mafft = ["mafft --maxiterate 1000 --genafpair --thread 10 --quiet", file, ">", file+'.mfa']
        if debug:
            print >>sys.stderr, ' '.join(mafft)
        os.system(' '.join(mafft))
    else:
        os.system("/bin/cp -f "+file+' '+file+'.mfa')
    #------------------------------------------------------
    psiblast = ['psiblast', '-db', target_db, '-out', file+'.psiblast.table', 
            '-evalue', evalue, '-outfmt', 
            '"7 qseqid sseqid qlen slen length pident qcovs evalue staxid ssciname scomname"', 
            '-num_threads 50 -num_iterations 10 -in_msa', file+'.mfa']
    if debug:
        print >>sys.stderr, ' '.join(psiblast)
    os.system(' '.join(psiblast))

    output = file+'.psiblast.table.parsed.tsv'
    parsePsiblastResultCMD(file+'.psiblast.table',identity,coverage,output)
    #return
    
    cmd = ["tail -n +2", output, "| cut -f 1 >", output+'.id']
    if debug:
        print ' '.join(cmd)
    os.system(' '.join(cmd))
    cmd = ["blastdbcmd -entry_batch", output+'.id', "-dbtype prot -db", target_db, "-out", output+'.fa']
    if debug:
        print ' '.join(cmd)
    #os.system(' '.join(cmd))
    
    #aDict = parsePsiblastResult(file+'.psiblast.table', identity, coverage)

    #final = file+'.psiblast.table.parsed.tsv'
    #final_fh = open(final, 'w')
    #outputPsiResult(aDict, fh=final_fh)
    #final_fh.close()

    ###--------multi-process------------------
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


