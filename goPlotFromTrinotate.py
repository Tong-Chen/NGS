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
    This is designed to summarize GO annotation results.
'''

import sys
import os
from json import dumps as json_dumps
from time import time, ctime,localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from math import log10

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
        metavar="FILEIN", help="The output of \
`extract_GO_assignments_from_Trinotate_xls.pl` \
with `include_ancestral_terms` turnned on. \
Normally `Trinotate/prefix.gene.GO`.")
    parser.add_option("-g", "--goLevel", dest="goLevel",
        metavar="FILEIN", default="/disk2/resource/goLevel.output", 
        help="go level file (default /disk2/resource/goLevel.output)")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_false", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def getGoLevel(goLevel):
    subrootD = {}
    sql = """
SELECT DISTINCT
  term.term_type, 
  term.acc, 
  term.name, 
  max(p.distance)
FROM
  term
  INNER JOIN graph_path AS p ON (p.term2_id=term.id)
  INNER JOIN term AS root ON (p.term1_id=root.id)
WHERE
  root.is_root=1
  AND term.is_obsolete=0
GROUP BY
  term.term_type, 
  term.acc, 
  term.name;
"""
    if not (goLevel or (os.path.exists("goLevel.output") and os.stat("goLevel.output").st_size>0)):
        file = 'goLevel.sql'
        file_fh = open(file, 'w')
        print >>file_fh, sql
        file_fh.close()
        #cmd = "mysql -A -hmysql-amigo.ebi.ac.uk -ugo_select \
#-pamigo -P4085 go_latest <goLevel.sql \
#| grep -P '\t2$' >goLevel.output"
        cmd = "mysql -A -hspitz.lbl.gov -ugo_select \
go_latest <goLevel.sql \
| grep -P '\t2$' >goLevel.output"
        print cmd
        #sys.exit(1)
        os.system(cmd)
        goLevel = "goLevel.output"
    #--------------------------------------------
    for line in open(goLevel):
        lineL = line.split("\t")
        subrootD[lineL[1]] = [lineL[0], lineL[2]]
    return subrootD
    
#----------------getGoLevel----------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    goLevel = options.goLevel
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    subrootD = getGoLevel(goLevel)
    goD = {}
    for line in fh:
        lineL = line.strip('\n').split('\t')
        goIdL = lineL[1].split(',')
        for goId in goIdL:
            if goId in subrootD:
                go_category, go_term = subrootD[goId]
                gokey = go_term
                if go_category not in goD:
                    goD[go_category] = {}
                if gokey not in goD[go_category]:
                    goD[go_category][gokey] = 1
                else:
                    goD[go_category][gokey] += 1
        #---------------------------------------------------------------
    #----------------------------------------------------------
    itemL = []
    go_categoryL = goD.keys()
    go_categoryL.sort()
    for go_category in go_categoryL:
        innerD = goD[go_category]
        for gokey, cnt in innerD.items():
            itemL.append([go_category, gokey, cnt])
        #--------------------------------------------------
    #--------------------------------------------------
    #print >>sys.stderr, itemL
    #itemL.sort(key=lambda x: x[2], reverse=True)
    #itemL = itemL[:top]
    itemL.sort(key=lambda x: (x[0], x[2]))
    count = []
    output = file + '.go_subparent.xls' 
    output_fh = open(output, 'w')
    print >>output_fh, "variable\tset\tvalue"
    for item in itemL:
        count.append(item[2])
        item[2] = str(item[2])
        print >>output_fh, '\t'.join(item)
    output_fh.close()

    maxCount = max(count)
    cntZero = int(log10(maxCount))
    y_label = ['0']
    for i in range(cntZero):
        y_label.append(str(10**(i+1)))
    y_label.append(str(maxCount))
    
    y_label = ','.join(y_label)

    levelL = ["\'"+item[1]+"\'" for item in itemL ] 
    cmd = "s-plot barPlot -f %s -d dodge -m TRUE -a set \
-R 90 -w 23 -u 18 -L \"%s\" -S 1 \
-y 'Number of Unigenes in each category' \
-v \"scale_y_log10(breaks=c(%s), labels=c(%s))\"" \
% (output, ','.join(levelL), y_label, y_label)
    #print cmd
    os.system(cmd)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
    if verbose:
        print >>sys.stderr,\
            "--Successful %s" % ctime()

if __name__ == '__main__':
    startTime = ctime()
    main()
    endTime = ctime()
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


