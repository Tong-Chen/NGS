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
    This is designed to do the GO enrichment in house.

The format of Gene Ontology file used in this program (could be generated using <generate_go_anno_for_goEnrichment.py>):

#--------Below is the file content (header line is needed)--------------------------
GO term    GO description  Gene_list   No of Genes under this GO term  Total annotated genes
GO:0045116  Function description	Bra002219,Bra006497	2	22567
GO:0004091	Function description    Bra000230,Bra000969,Bra004627,Bra012745,Bra018016,Bra021530,Bra028144	7	22567
GO:0004420	Function description    Bra002053,Bra008261,Bra014102,Bra015739	4	22567
GO:0004830	Function description    Bra001114,Bra028092,Bra029108,Bra034277,Bra040117,Bra040132	6	22567
#--------Above is the file content--------------------------

'''

import sys
import os
from fisher import pvalue
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from statsmodels.stats.multitest import multipletests

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
    parser.add_option("-i", "--input-file", dest="go",
        metavar="go", help="The gene ontology file")
    parser.add_option("-g", "--gene", dest="gene",
        metavar="gene", help="One column gene list file or multiple columns only the first column will be used.")
    parser.add_option("-p", "--pvalue", dest="pvalue",
        default=0.05, metavar="pvalue", help="pvalue for enriched terms. Default 0.05.")
    parser.add_option("-o", "--output", dest="output",
        help="Specify output file name if multiple-test correction is needed. Default STDOUT. If there are 2 columns in gene list file, the second column will be used as output filename. Given string will be used as prefix.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.go != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    go_file = options.go
    gene_file = options.gene
    output = options.output
    pvalue_thresh = float(options.pvalue)
    verbose = options.verbose
    debug = options.debug
    #-----------------------------------
    annoGeneD = set([j for i in open(go_file) \
        for j in i.split('\t')[2].split(',')])
    #-----------------------------------------------
    geneGrpD = {}
    for line in open(gene_file):
        lineL = line.strip().split('\t')
        len_lineL = len(lineL)
        if len_lineL == 2:
            gene, grp = lineL
            grp = grp.replace(' ','_')
        elif len_lineL == 1:
            gene = lineL[0]
            grp = 1
        if grp not in geneGrpD:
            geneGrpD[grp] = set([])
        geneGrpD[grp].add(gene)
    #----------------------------------------------
    #geneD = set([i.split("\t")[0].strip() for i in open(gene_file)])
    for grp, geneD in geneGrpD.items():
        if grp == 1:
            if output:
                fh_out = open(output, 'w')
            else:
                fh_out = sys.stdout
        else:
            if output:
                fh_out = open(output+'.'+grp, 'w')
            else:
                fh_out = open(gene_file+'.'+grp,'w')
        #---------------------------------------------
        geneD = geneD.intersection(annoGeneD)
        geneD_len = len(geneD)
        #--------------------------------
        annoL = []
        header = 1
        tmpL = []
        pL = []
        headerL = []
        for line in open(go_file):
            lineL = line.strip().split('\t')
            if header:
                headerL = lineL[:]
                #print >>fh_out, "%s\t%s\tTargetGene\tTargetCount\tTargetTotal\tPvalue\tOdds ratio\tFDR" \
                #    % (lineL[0], lineL[1])
                header -= 1
                continue
            go_gene = lineL[2].split(',')
            anno_gene = [gene for gene in go_gene if gene in geneD]
            anno_gene = set(anno_gene)
            annoCount = len(anno_gene)
            if annoCount > geneD_len:
                print >>sys.stderr, go_gene
                print >>sys.stderr, anno_gene
                sys.exit(1)
            if annoCount > 3:
                termCount = int(lineL[3])
                totalCount = int(lineL[4])
                annoCount = len(anno_gene)
                #print >>sys.stderr, termCount, totalCount, annoCount
                p = pvalue(annoCount, termCount-annoCount,
                        geneD_len-annoCount,
                        totalCount-geneD_len+annoCount-termCount)
                p = p.right_tail
                fracT = annoCount * 1.0 / geneD_len / termCount * totalCount
                if fracT > 1 and p <= pvalue_thresh*2:
                    fracT = "%.2f" % fracT
                    tmpL.append([lineL[0], lineL[1], ', '.join(anno_gene), str(annoCount), 
                        str(geneD_len), p, fracT,str(p)])
                    pL.append(p) 
        #-------------END reading file----------
        #tmpL.sort(key=lambda x: float(x[5]))
        try:
            p_adjL = multipletests(pL, method="fdr_bh")[1]
        except ZeroDivisionError:
            print >>sys.stderr, pL
            print >>sys.stderr, tmpL
            continue
        #---------------------------------------------
        for eachtmpL, p_adj in zip(tmpL, p_adjL):
            eachtmpL[-3] = format(eachtmpL[-3],'0.2E')
            eachtmpL[-1] = format(p_adj,'0.2E')
        #------------------------------------------------
        tmpL.sort(key=lambda x: float(x[5]))

        for i in tmpL:
            if float(i[-1]) > pvalue_thresh:
                break
            if headerL:
                print >>fh_out, "%s\t%s\tTargetGene\tTargetCount\tTargetTotal\tPvalue\tOdds ratio\tFDR" \
                    % (headerL[0], headerL[1])
                headerL = []
            print >>fh_out, '\t'.join(i)

        #print >>fh_out, '\n'.join(['\t'.join(i) for i in tmpL])
        #----------------------
        if output or grp != 1:
            fh_out.close()
            #os.system("multipleTest.sh -f %s" % output)
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


