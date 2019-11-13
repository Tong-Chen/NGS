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
    This is designed to integrate the results of <kmeans.test.sh>,
    and <runGOseq.pl>.


go_assignments (first file):
    ENSOARG00000014536      GO:0030529, GO:0006412, GO:0005840,
    ENSOARG00000007374      GO:0006954, GO:0008284, GO:0001568,
    ENSOARG00000014537      GO:0005515, GO:0005515, GO:0003676,
    ENSOARG00000007375      GO:0005515, GO:0002376,

go_assignments (second file, only first two columns will be used):
    ENSOARG00000014536  ATGC    123
    ENSOARG00000007374  OCT4    786
    ENSOARG00000014537  WOX     3237237
    ENSOARG00000007375  IDD     2372

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
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="factor_labeling given to runGOSeq.pl.")
    parser.add_option("-D", "--dir", dest="dir",
            default='.', 
        metavar="DIRECTORY", help="The directory containing all data.\
normally \
<quantification/trinity_quantification.DE/$(prefix).matrix.trend.16.kmeans>.")
    parser.add_option("-o", "--out_prefix", dest="out_prefix",
        help="Default the out_prefix would be the innest dir \
given to <-D>.")
    parser.add_option("-g", "--go_assignments", dest="go_assign",
        help="GO_assignments file given to run_GOseq.pl. Zero or \
one or two files separated by ',' is accepted. \
If two files are given, the second file may contain at least \
two columns. The first column should be gene names of the first \
file. The second column may be gene names in other format. \
Other columns will be ignored.")
    parser.add_option("-p", "--pvalue", dest="pvalue",
        default=0.05, help="FDR used for selecting DE genes")
    parser.add_option("-f", "--log2fc", dest="log2fc",
        default=2, help="log2 fold change used for selecting DE genes")
    parser.add_option("-e", "--expr", dest="expr",
        help="<quantification/$(prefix).isoform.expr.TMM.fpkm.matrix>")
    parser.add_option("-a", "--anno", dest="anno",
        help="<Trinotate/Trinotate_annotation_report.xls>")
    parser.add_option("-I", "--annoIndex", dest="index",
        default=2, help="For <isoforms>, index should be 2 \
representing indexing annotation by the second column ids. \
For <genes>, index should be 1 meaning indexing annotation \
by the first column.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def readEnriched(*enriched):
    header = 1
    for label, file, file_out, goAD, goNameD, geneL in enriched:
        label += "-UP.enriched"
        out_fh = open(file_out, 'w')
        header = 1
        for line in open(file):
            lineL = line.strip().split('\t')
            if len(lineL) < 9:
                print >>sys.stderr, "%s NO ENRICHEMNT" % file
                break
            goNum = lineL[0]
            #print >>sys.stderr, file
            #print >>sys.stderr, lineL
            newLineL = [goNum,lineL[8], lineL[1], lineL[3], lineL[7]]

            if header:
                newLineL.insert(0, "Sample")
                newLineL.append("Gene_name")
                newLineL.append("Gene_alias_name")
                header -= 1
            else:
                '''
                goAD = {
                    # In this script,  only file1 is used, no file2 in goAD
                    GO1:{file1:[geneA, geneB], file2:[geneA, geneB]}, 
                    GO2:{file1:[geneC, geneD], file2:[geneE, geneF]}, 
                }'''
                goD = goAD.get(goNum, "")
                if goD:
                    #valueL = [','.join(i) for i in goD.values() \
                    #    for j in i if j in geneL]
                    valueL = []
                    for i in goD.values():
                        tmpL = []
                        for j in i:
                            if j in geneL:
                                tmpL.append(j)
                        valueL.append(','.join(tmpL))
                    #----------------------------------
                    tmpL = valueL[-1].split(',')
                    valueL.append(','.join(\
                        [goNameD.get(i, " ") for i in tmpL]))

                    newLineL.extend(valueL)
                newLineL.insert(0, label)
            print >>out_fh, '\t'.join(newLineL)
        out_fh.close()
    #--------------------------------------
#--------------------------------------

def readDepleted(*depleted):
    header = 1
    for label, file, file_out, goAD in depleted:
        label += "-UP.depleted"
        out_fh = open(file_out, 'w')
        header = 1
        #print >>sys.stderr, file
        for line in open(file):
            lineL = line.strip().split('\t')
            #print >>sys.stderr, lineL
            goNum = lineL[0]
            if len(lineL) < 10:
                break
            newLineL = [goNum, lineL[9], lineL[2], lineL[3], lineL[8], lineL[9]]
            if header:
                newLineL.insert(0,"Sample")
                header -= 1
            else:
                goD = goAD.get(goNum, "")
                if goD:
                    valueL = [','.join(i) for i in goD.values()]
                    newLineL.extend(valueL)
                newLineL.insert(0,label)
            print >>out_fh, '\t'.join(newLineL)
        out_fh.close()
    #--------------------------------------
#--------------------------------------

def getTopTerm(type, prefix, top, *file):
    maxLen = 70
    file_out_n = prefix+".GOseq."+type+'.xls'
    file_out = open(file_out_n, 'w')
    print >>file_out, "Sample\tGO\tTerm\tneg_log10pvalue\tCount\tFDR\tGene"
    for single in file:
        count = 0
        for line in open(single):
            if line.find("\t"+type) != -1:
                lineL = line.split('\t')[:7]
                lineL[0] = lineL[0][:maxLen]
                print >>file_out, '\t'.join(lineL)
                count += 1
            if count >= top:
                break
    #---------------------------------------
    file_out.close()
    if count:
        height = count / 3
        if height < 15:
            height = 15
        elif height < 25:
            height = 25
        cmd = ['s-plot scatterplotDoubleVariable -f', file_out_n, 
            '-o Sample -v Term -c neg_log10pvalue -s Count -w 25 -a', 
            str(height), '-E pdf -R 30 -H 1 -V 1 -l neg_log10pvalue']
        os.system(' '.join(cmd))
        convert = ['convert -density 150 -quality 90',
                file_out_n+'.scatterplot.dv.pdf',
                file_out_n+'.scatterplot.dv.png']
        os.system(' '.join(convert))
#------------------------------------

def readAnno(anno, index):
    header = 1
    annoD = {}
    for line in open(anno):
        line = line.strip()
        if header:
            head = line
            header -= 1
        else:
            lineL = line.split('\t', 3)
            annoD[lineL[index]] = line
    #--------------------------
    return annoD, head
#-----------------------------
def readExpr(expr):
    header = 1
    exprD = {}
    for line in open(expr):
        line = line.rstrip()
        if header:
            head = 'expr_TMM_FPKM'+line
            header -= 1
        else:
            lineL = line.split('\t', 2)
            exprD[lineL[0]] = line
    #----------------
    return exprD, head
#---------------------------
def annoCluster(geneL, annoD, annoH, exprD, exprH, file_out):
    file_fh = open(file_out, 'w')
    #print >>file_fh, '\t'.join(["ID", exprH, annoH])
    print >>file_fh, '\t'.join([exprH, annoH])
    for key in geneL:
        #print >>file_fh, '\t'.join([exprD[key], annoD[key]])
        print >>file_fh, '\t'.join([exprD.get(key, ""), annoD.get(key, '')])
    file_fh.close()
#_----------------------------------------

def annoSubset(annoD, annoH, *subset):
    for file in subset:
        file_out =  file.replace('subset', 'xls')
        file_fh = open(file_out, 'w')
        header = 1
        for line in open(file):
            line = line.strip()
            if header:
                print >>file_fh, '\t'.join([line, annoH])
                header -= 1
            else:
                key = line.split('\t', 1)[0]
                print >>file_fh, '\t'.join([line, annoD[key]])
        file_fh.close()
#_----------------------------------------

#def readGo_assign(go_assign):
#    '''
#    go_assign = "file1,file2,file3"
#
#    input file:
#        geneA   GO1,GO2,GO3
#        geneB   GO1,GO2,GO3,GO4
#    
#    goAD = {
#        GO1:{file1:[geneA, geneB], file2:[geneA, geneB]}, 
#        GO2:{file1:[geneC, geneD], file2:[geneE, geneF]}, 
#    }
#    '''
#    goAD = {}
#    for file in go_assign.split(','):
#        name = os.path.split(file)[1]
#        assert name not in goAD, "Duplicate %s" % name
#        #goAD[name] = {}
#        for line in open(file):
#            gene, go = line.strip().split('\t')
#            goL = go.split(',')
#            for go in goL:
#                if go not in goAD:
#                    goAD[go] ={}
#                if name not in goAD[go]:
#                    goAD[go][name] = []
#                goAD[go][name].append(gene)
#            #-----------------------------------
#        #---------------------------------------
#    return goAD
##-------------END goAD--------------------


def readGo_assign(go_assign):
    '''
    go_assign = "file1,file2"

    input file1:
        geneA   GO1,GO2,GO3
        geneB   GO1,GO2,GO3,GO4
    
    input file2 (only first two columns will be used):
        geneA   symbolA entrezA
        geneB   symbolB entrezB
        geneC   symbolC entrezC

    goAD = {
        GO1:{file1:[geneA, geneB], file2:[geneA, geneB]}, 
        GO2:{file1:[geneC, geneD], file2:[geneE, geneF]}, 
    }

    goNameD = {
        geneA: [symbolA, entrezA], 
        geneB: [symbolB, entrezB], 
        geneC: [symbolC, entrezC]
    }

    '''
    goAD = {}
    fileL = go_assign.split(',')
    file = fileL[0]
    name = os.path.split(file)[1]
    assert name not in goAD, "Duplicate %s" % name
    #goAD[name] = {}
    for line in open(file):
        gene, go = line.strip().split('\t')
        goL = go.split(',')
        for go in goL:
            if go not in goAD:
                goAD[go] ={}
            if name not in goAD[go]:
                goAD[go][name] = []
            goAD[go][name].append(gene)
        #-----------------------------------
    #---------------------------------------
    goNameD = {}
    if (len(fileL) > 1):
        file = fileL[1]
        for line in open(file):
            lineL = line.strip('\n').split('\t')
            key = lineL[0]
            goNameD[key] = lineL[1]
    #-------------------------------------
    return goAD, goNameD
#-------------END goAD--------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    factor_labeling = options.filein
    dir = options.dir + '/'
    out_prefix = options.out_prefix
    if not out_prefix:
        if dir == './':
            out_prfix = ''
        else:
            out_prefix = dir + dir.rstrip('/').split('/')[-1] + '.'
    verbose = options.verbose
    pvalue = options.pvalue
    log2fc = options.log2fc
    global debug
    debug = options.debug
    anno = options.anno
    index = int(options.index) - 1
    if anno:
        annoD, annoH = readAnno(anno, index)
    #print >>sys.stderr, annoD['TR6117|c0_g1_i1']
    expr = options.expr
    if expr:
        exprD, exprH = readExpr(expr)
    #print >>sys.stderr, exprD['TR6117|c0_g1_i1']
    go_assign = options.go_assign
    goAD = {}
    if go_assign:
        goAD, goNameD = readGo_assign(go_assign)
    #-----------------------------------
    if factor_labeling == '-':
        fh = sys.stdin
    else:
        fh = open(factor_labeling)
    #--------------------------------
    factorL = []
    factorD = {}
    for line in fh:
        gene, factor = line.strip().split('\t')
        if factor not in factorD:
            factorD[factor] = []
            factorL.append(factor)
        factorD[factor].append(gene)
    #----------------------------------------------
    for factor in factorL:
        geneL = factorD[factor]
        output = out_prefix + factor + '.anno.xls'
        if anno and expr:
            annoCluster(geneL, annoD, annoH, exprD, exprH, output)
        
        enriched = dir + factor + '.GOseq.enriched'
        enriched_out = out_prefix+factor+'.GOseq.enriched.xls'
        depleted = dir + factor + '.GOseq.depleted'
        depleted_out = out_prefix+factor+'.GOseq.depleted.xls'

        readEnriched([factor, enriched, enriched_out, goAD, goNameD, geneL])
        #readDepleted([factor, depleted, depleted_out, goAD])
        
        getTopTerm('BP', out_prefix+factor, 20, enriched_out)
        getTopTerm('MF', out_prefix+factor, 20, enriched_out)
        getTopTerm('CC', out_prefix+factor, 20, enriched_out)
        #getTopTerm('BP', out_prefix+factor, 50, enriched_out, depleted_out)
        #getTopTerm('MF', out_prefix+factor, 50, enriched_out, depleted_out)
        #getTopTerm('CC', out_prefix+factor, 50, enriched_out, depleted_out)

        #getTopTerm('MF', prefix, 50, subsetA_enriched+'.xls',
        #    subsetA_depleted+'.xls', subsetB_enriched+'.xls',
        #    subsetB_depleted+'.xls')

        #getTopTerm('CC', prefix, 50, subsetA_enriched+'.xls',
        #    subsetA_depleted+'.xls', subsetB_enriched+'.xls',
        #    subsetB_depleted+'.xls')
    #-------------END reading factor_labeling----------
    #----close factor_labeling handle for files-----
    if factor_labeling != '-':
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
            "--Successful %s" % strftime(timeformat, localtime())

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()
    ###---------procompare_pair the program---------
    #import procompare_pair
    #procompare_pair_output = sys.argv[0]+".prof.txt")
    #procompare_pair.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(procompare_pair_output)
    #p.sort_stats("time").print_stats()
    ###---------procompare_pair the program---------


