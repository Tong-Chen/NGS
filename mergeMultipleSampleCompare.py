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
   This is designed to merge multipleSampleCompare result. 

Parameters:
    
    -i: should be a string in format like 
        <bin50Cnt.merge/SCX.CTCTCT.bin50Cnt.merge.all.DE>.
        Actually this is the summary file output by <multipleSampleCompare.py>.
        Every word should be changed accordingly except <CTCTCT> which will be
        substituted by strings given to <-l>.

        The <bin50Cnt.merge/SCX.CTCTCT.bin50Cnt.merge.all.DE>, 
        <bin50Cnt.merge/SCX.CTCTCT.bin50Cnt.merge.all.DE.norm>, 
        <bin50Cnt.merge/SCX.CTCTCT.bin50Cnt.merge.all.DE.count.xls> files would be used.

        Besides the differentially expressed or modified files will be extracted 
        based on the second column of <bin50Cnt.merge/SCX.CTCTCT.bin50Cnt.merge.all.DE> 
        and will be merged together.

        Output files include <bin50Cnt.merge/SCX.bin50Cnt.merge.all.DE>, 
        <bin50Cnt.merge/SCX.bin50Cnt.merge.all.DE.norm>, 
        <bin50Cnt.merge/SCX.bin50Cnt.merge.all.DE.count.xls>.

    -l: list of separated files. When designing this program, we separated whole
        genome data into each chromosome. So here, the list of chormosomes should
        be supplied in format like 'chr1 chr2 chr3 .. chrY' (blank spearated).
        These strings will be used to substitute <CTCTCT> in <-i> to get reak filename.

    -s: strings as given to <-i>

A series of files may be needed are:

SCX.chr4.bin50Cnt.merge.all.DE
''
chr4_3803596    OA_SF_280._higherThan_.RA_SF_280
chr4_3803595    OA_SF_280._higherThan_.RA_SF_280
chr4_3803594    OA_SF_280._higherThan_.RA_SF_280
chr4_981852     OA_PM._lowerThan_.OA_SF_350
chr4_982984     OA_PM._lowerThan_.OA_SF_350
chr4_982985     OA_PM._lowerThan_.OA_SF_350
''

SCX.chr4.bin50Cnt.merge.RA_SF_280._lowerThan_.RA_SF_450.anno.xls
SCX.chr4.bin50Cnt.merge.RA_SF_280._higherThan_.RA_SF_450.anno.xls

SCX.chr4.bin50Cnt.merge.RA_SF_350._vs_.RA_SF_450.anno.xls

''
This needs compare_pair file
''

SCX.chr4.bin50Cnt.merge.RA_SF_350._lowerThan_.RA_SF_450.anno.xls
SCX.chr4.bin50Cnt.merge.RA_SF_350._higherThan_.RA_SF_450.anno.xls
SCX.chr4.bin50Cnt.merge.all.DE.count.xls
''
Samp OA_PM OA_SF OA_SF_280 OA_SF_350 OA_SF_450 RA_PM RA_SF RA_SF_280 RA_SF_350 RA_SF_450
OA_PM NA NA 45 175 NA NA NA NA NA NA
OA_SF NA NA NA NA NA NA NA NA NA NA
OA_SF_280 NA NA NA NA NA NA NA NA NA NA
OA_SF_350 1 NA NA NA NA NA NA NA NA NA
OA_SF_450 NA NA NA NA NA NA NA NA NA NA
RA_PM NA NA NA NA NA NA 24 NA NA NA
''
SCX.chr4.bin50Cnt.merge.all.DE.count.xls.heatmapS.pdf
SCX.chr4.bin50Cnt.merge.all.DE.norm
    
'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
import gzip


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
        metavar="FILEIN", help="A string as format specified above.")
    parser.add_option("-l", "--list", dest="list",
        help="A list of strings as described above.")
    parser.add_option("-c", "--compare-pair", dest="compare_pair",
        help="<compare_pair> file given to multipleSampleCompare.py.")
    parser.add_option("-s", "--other-files-merge", dest="file_merge",
        help="<bin50Cnt.merge/SCX.CTCTCT.bin50Cnt.merge.norm.xls.gz.std.top300000.xls.gz>")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file_str  = options.filein
    subL  = options.list.split()
    file_merge = options.file_merge
    verbose = options.verbose
    compare_pair = options.compare_pair
    compare_pairL = [line.strip().replace('\t', '._vs_.') for line in open(compare_pair)]
    global debug
    debug = options.debug
    #-----------------------------------
    annoD = {}
    all_de = file_str.replace('CTCTCT.', '')
    all_de_fh = open(all_de, 'w')
    all_de_norm = all_de+'.norm'
    all_de_norm_fh = open(all_de_norm, 'w')
    de_region_countD = {}
    all_de_count = all_de + '.count.xls'
    count_header = 1
    for sub in subL:
        #file_str would be <bin50Cnt.merge/SCX.CTCTCT.bin50Cnt.merge.all.DE>
        #sub_file_de would be <bin50Cnt.merge/SCX.chr4.bin50Cnt.merge.all.DE>
        sub_file_de = file_str.replace('CTCTCT', sub)
        #sub_file_list = [i.split()[1] for i in open(sub_file_de)]
        #sub_file_list = list(set(sub_file_list))
        sub_file_list = set()
        for de_line in open(sub_file_de):
            print >>all_de_fh, de_line, 
            sub_file_list.add(de_line.split()[1])
        #de norm
        sub_file_de_norm = sub_file_de + '.norm'
        sub_file_de_norm_header = 1
        for line in open(sub_file_de_norm):
            if sub_file_de_norm_header:
                if count_header:
                    print >>all_de_norm_fh, line, 
                    count_header -= 1
                sub_file_de_norm_header -= 1
                continue
            #----------------------------------
            print >>all_de_norm_fh, line, 

        #DE count
        sub_file_de_count_header = 1
        sub_file_de_count = sub_file_de + '.count.xls'
        for de_count_line in open(sub_file_de_count):
            de_count_lineL = de_count_line.split()
            if sub_file_de_count_header:
                de_count_headerL = de_count_lineL
                sub_file_de_count_header -= 1
                continue
            key = de_count_lineL[0]
            if key not in de_region_countD:
                de_region_countD[key] = {}
            for tmpkey, tmpvalue in zip(de_count_headerL[1:], de_count_lineL[1:]):
                if tmpvalue != 'NA':
                    tmpvalue = int(tmpvalue)
                if tmpkey not in de_region_countD[key]:
                    de_region_countD[key][tmpkey] = tmpvalue
                else:
                    if de_region_countD[key][tmpkey] == 'NA' and tmpvalue != 'NA':
                        de_region_countD[key][tmpkey] = tmpvalue
                    else:
                        #print >>sys.stderr, de_region_countD
                        #print >>sys.stderr, de_region_countD[key]
                        #print >>sys.stderr, de_region_countD[key][tmpkey]
                        #print >>sys.stderr, tmpvalue
                        if tmpvalue != 'NA':
                            de_region_countD[key][tmpkey] += tmpvalue
            #-----------------------------------------
        #--------DE count---------------------------------

        #sub_anno_prefix would be <bin50Cnt.merge/SCX.chr4.bin50Cnt.merge.>
        #Subfiles would be generated by combining this prefix and the second 
        #column in this file
        sub_anno_prefix = file_str.replace('all.DE', '')
        #anno_prefix = file_str.replace('CTCTCT.', '').replace('all.DE', '')
        if debug:
            print >>sys.stderr, sub, sub_file_list
        #for compare_pair in compare_pairL:
        #    sub_vs_input = sub_anno_prefix.replace('CTCTCT',sub)+compare_pair+'.anno.xls'
        #    sub_vs_input = sub_anno_prefix.replace('CTCTCT',su)+compare_pair+'.anno.xls'
        #----------------------------------------------
        sub_file_list = list(sub_file_list)
        sub_file_list.extend(compare_pairL)
        for sub_file in sub_file_list:
            sub_anno_input  = sub_anno_prefix.replace('CTCTCT',sub)+sub_file+'.anno.xls'
            sub_anno_output = sub_anno_prefix.replace('CTCTCT.','')+sub_file+'.anno.xls'
            if debug:
                print >>sys.stderr, sub_anno_input
                print >>sys.stderr, sub_anno_output
            header = 1
            for line in open(sub_anno_input):
                if header:
                    header -= 1
                    if sub_anno_output not in annoD:
                        annoD[sub_anno_output] = open(sub_anno_output, 'w')
                        print >>annoD[sub_anno_output], line,
                    continue
                #---------------------------------
                print >>annoD[sub_anno_output], line,
            #-------------------------------------
        #------------------------------------------
    all_de_fh.close()
    all_de_norm_fh.close()
    for value_fh in annoD.values():
        value_fh.close()
    ###--------multi-process------------------
    file_merge_output = file_merge.replace('CTCTCT.', '').replace('.gz', '')
    file_merge_output_fh = open(file_merge_output, 'w')
    file_merge_output_header = 1
    for sub in subL:
        sub_file_merge = file_merge.replace('CTCTCT', sub)
        header = 1
        for line in gzip.open(sub_file_merge, 'rb'):
            if header:
                if file_merge_output_header:
                    print >>file_merge_output_fh, line
                    file_merge_output_header -= 1
                header -= 1
                continue
            print >>file_merge_output_fh, line, 
        #---------------------------------
    file_merge_output_fh.close()
    
    #cmdL = ['s-plot pca -f', file_merge_output]
    #os.system(' '.join(cmdL))

    #de_region_countD = {}
    #all_de_count = all_de + '.count.xls'
    #de_count_headerL
    all_de_count_fh = open(all_de_count, 'w')
    print >>all_de_count_fh, '\t'.join(de_count_headerL)
    for samp in de_count_headerL[1:]:
        tmpL = [str(de_region_countD[samp][i]) for i in de_count_headerL[1:]]
        print >>all_de_count_fh, "{}\t{}".format(samp, '\t'.join(tmpL))
    all_de_count_fh.close()
    all_de_count_cmd = ['s-plot heatmapS -f', all_de_count, "-A 45 -T 2 -l top -I 'DE genes or regions count'", "-b TRUE -Y white"]
    os.system(' '.join(all_de_count_cmd))
    #-----------------------------------------

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


