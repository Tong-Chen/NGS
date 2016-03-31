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
    This is used to statistics tophat/bowtie2 mapping quality.
'''

import sys
import re
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

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
    parser.add_option("-i", "--input-files", dest="filein",
        metavar="FILEIN", help="The output of tophat/bowtie2 \
align_summary.txt. \
Multiple files can be separated with \
<,> or space(< >) or both, like <'file1,file2,file3'>, \
<'file1 file2 file3'> or <'file1' 'file2', 'file3'>." )
    parser.add_option("-l", "--file-label", dest="file_label",
        metavar="FILE_LABEL", help="The labels should have the same \
format and order as filenames given to <-i>.")
    parser.add_option("-S", "--seq-type", dest="seq_type",
        default='PE', help="PE or SE")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    fileL = re.split('[, ]*', options.filein.strip())
    labelL = re.split('[, ]*', options.file_label.strip())
    lenfileL = len(fileL)
    lenlabeL = len(labelL)
    assert lenfileL == lenlabeL
    seq_type = options.seq_type
    debug = options.debug
    #-----------------------------------
    #print "# Mapping quality"
    #print 
    #print "## Mapping quality"
    #print 

    #print "Table: Statistics of mapping status for all samples\n"

    outputL = []
    dict = {}
    #-----------------------------------
    for i in range(lenfileL):
        file = fileL[i]
        fh = open(file)
        label = labelL[i]
        if seq_type == 'PE':
            line = fh.readline() #Left reads
            left_tag = int(fh.readline().split(':')[1].strip())
            #print left_tag
            left_map = int(fh.readline().split(':')[1].strip().split()[0])
            #print left_map
            left_multiple = int(fh.readline().split(':')[1].strip().split()[0])
            #print left_multiple

            line = fh.readline() #Right reads
            right_tag = int(fh.readline().split(':')[1].strip())
            assert left_tag == right_tag, ""
            #print right_tag
            right_map = int(fh.readline().split(':')[1].strip().split()[0])
            #print right_map
            right_multiple = int(fh.readline().split(':')[1].strip().split()[0])
            #print right_multiple
            
            line = fh.readline() #overall read mapping rate
            line = fh.readline() # blank line
            aligned_pair = int(fh.readline().split(':')[1].strip())
            #print aligned_pair
            line = fh.readline()
            discordant_align_pair = int(fh.readline().strip().split()[0])
            #print discordant_align_pair

            totaltag = left_tag + right_tag
            total_map_tag = left_map + right_map
            total_map_tag_per = "%.1f" % (total_map_tag * 100.0 / totaltag)
            total_unique_tag = total_map_tag - left_multiple - right_multiple
            total_unique_tag_per = "%.1f" % (total_unique_tag * 100.0 / totaltag)
            
            print "%s.map=%d" % (label, total_unique_tag)
#            aligned_pair_per = "%.1f" % (aligned_pair * 100.0 / left_tag)
#            aligned_concordant_pair = aligned_pair - discordant_align_pair
#            aligned_concordant_pair_per = "%.1f" % \
#                (aligned_concordant_pair * 100.0 / left_tag)
#            sampL_out = [label, str(totaltag),
#                str(total_map_tag)+' ('+str(total_map_tag_per)+')', 
#                str(total_unique_tag)+' ('+str(total_unique_tag_per)+')', 
#                str(aligned_pair)+' ('+str(aligned_pair_per)+')', 
#                str(aligned_concordant_pair)+' ('+str(aligned_concordant_pair_per)+')'
#            ]
#            outputL.append(sampL_out[:])
        else: #SE
            line = fh.readline() #Reads reads
            line = fh.readline() #Input
            line = fh.readline() #mapped
            total_unique_tag = line.split(':')[1].strip().split()[0]
            print "%s.map=%s" % (label, total_unique_tag)
    #-------------END reading file----------
    #print '\n'.join(['\t'.join(i) for i in outputL])
    #transform(outputL)

    #col = len(outputL)
    #row = len(outputL[0])
    
    #newOutPutL = []

    #for i in range(row):
    #    newOutPutL.append([outputL[j][i] for j in range(col)])
    
    #print '\n'.join(['\t'.join(i) for i in newOutPutL])
    #print '|'.join(headerL)
    #print '|'.join(slashL)
    #print '\n'.join(['|'.join(i) for i in outputL[1:]])
    #print

    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------

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


