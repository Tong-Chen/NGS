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

The format of Gene Ontology file used 

## Optional annotation format--(long, at least 3 columns)--------
## If there is a forth column, it will be treated as different database annotation, such
## as GO BP,  GO CC, KEGG , etc.
GO term GO Description  Gene
GO:0045116  Function description    Bra002219
GO:0045116  Function description    Bra006497
GO:0004091  Function description    Bra000230
GO:0004091  Function description    Bra000969
GO:0004091  Function description    Bra004627
#--------------------------------------

TEST WELL
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
    parser.add_option("-i", "--go-anno-file", dest="go",
        metavar="go", help="The gene ontology file containing 3 columns. The first two columns will be used together and can be any annotation or any string when empty. The third column must be gene name with each one line.")
    parser.add_option("-t", "--annoType", dest="annoType",
        default="Anno", help="Any string like BP, MF, CC or KEGG, Reactome represents what type of annotation used.")
    parser.add_option("-g", "--gene", dest="gene",
        metavar="gene", help="One column gene list file or two columns with additional group column for each gene.")
    parser.add_option("-p", "--pvalue", dest="pvalue",
        default=0.05, metavar="pvalue", help="pvalue for enriched terms. Default 0.05.")
    parser.add_option("-o", "--output", dest="output",
        help="Prefix for output files. If there are 2 columns in gene list file, the second column will be used as a tag for output filename. If there are annotations from different databases (4 columsn in annotation file), db names will be used as a tag for output filename also.")
    parser.add_option("-P", "--toptermsPlot", dest="topPlot",
        default = 0, type="int", 
        help="Specify number of top items for enrichment plot. Default 0. Accept integers like 20.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.go != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def readLongIn(file, label="Anno"):
    gene_col = 2
    go_id_col = 0
    go_desp_col = 1
    header = 1
    goAnnoD = {}
    geneSet = set([])
    allAnnoD = {}
    for line in open(file):
        if header:
            header -= 1
            continue
        lineL = line.strip().split('\t')
        gene  = lineL[gene_col]
        go_id = lineL[go_id_col]
        go_desp = lineL[go_desp_col]
        if (len(lineL)==4):
            label = lineL[3]
        if label not in allAnnoD:
            allAnnoD[label] = {}
            allAnnoD[label]['goAnnoD'] = {}
            allAnnoD[label]['geneSet'] = set([])
        if id == '--' and desp == '--':
            continue
        #geneSet.add(gene)
        allAnnoD[label]['geneSet'].add(gene)
        key = go_id+'\t'+go_desp
        if key not in allAnnoD[label]['goAnnoD']:
            allAnnoD[label]['goAnnoD'][key] = set([gene])
        else:
            allAnnoD[label]['goAnnoD'][key].add(gene)
    #-------------END reading file----------
    #totalGene = len(geneSet)
    #return goAnnoD, totalGene, geneSet
    return allAnnoD
#---------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    go_file = options.go
    gene_file = options.gene
    output = options.output
    pvalue_thresh = float(options.pvalue)
    verbose = options.verbose
    debug = options.debug
    topPlot = options.topPlot
    annoType = options.annoType
    #-----------------------------------
    #annoGeneD = set([j for i in open(go_file) \
    #    for j in i.split('\t')[2].split(',')])
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
    #goAnnoD, total_anno_gene, total_anno_gene_set = readLongIn(go_file, label=annoType)
    allAnnoD = readLongIn(go_file, label=annoType)
    
    for label, valueD in allAnnoD.items():
        goAnnoD = valueD['goAnnoD']
        total_anno_gene_set = valueD['geneSet']
        total_anno_gene = len(total_anno_gene_set)
        annoType = label
        out_file = ''
        for grp, geneD in geneGrpD.items():
            if grp == 1:
                if output:
                    out_file = output+'.'+label+'.xls'
                    fh_out = open(output, 'w')
                else:
                    fh_out = sys.stdout
            else:
                if output:
                    out_file = output+'.'+grp+'.'+label+'.xls'
                    fh_out = open(out_file, 'w')
                else:
                    out_file = gene_file+'.'+grp+'.'+label+'.xls'
                    fh_out = open(out_file,'w')
            #---------------------------------------------
            geneD = geneD.intersection(total_anno_gene_set)
            geneD_len = len(geneD)
            #--------------------------------
            annoL = []
            header = 1
            tmpL = []
            pL = []
            headerL = []
            for go_id, go_gene in goAnnoD.items():
                #go_desp, go_gene = itemL
                #anno_gene = [gene for gene in go_gene if gene in geneD]
                #anno_gene = set(anno_gene)
                anno_gene = go_gene.intersection(geneD)
                annoCount = len(anno_gene)
                if annoCount > geneD_len:
                    print >>sys.stderr, go_gene
                    print >>sys.stderr, anno_gene
                    sys.exit(1)
                if annoCount > 3:
                    termCount = len(go_gene)
                    p = pvalue(annoCount, termCount-annoCount,
                            geneD_len-annoCount,
                            total_anno_gene-geneD_len+annoCount-termCount)
                    #p = p.two_tail
                    p = p.right_tail
                    #print >>sys.stderr, [go_desp, annoCount-1, termCount-annoCount, geneD_len-annoCount, total_anno_gene-geneD_len+annoCount-termCount, p]
                    #print >>sys.stderr, [go_desp, annoCount, termCount, geneD_len, total_anno_gene, p]
                    fracT = annoCount * 1.0 / geneD_len / termCount * total_anno_gene
                    if fracT > 1 and p <= pvalue_thresh*2:
                        fracT = "%.2f" % fracT
                        tmpL.append([go_id, ', '.join(anno_gene), str(annoCount), 
                            str(geneD_len), str(termCount), str(total_anno_gene), fracT,p,str(p)])
                        pL.append(p) 
            #-------------END reading file----------
            try:
                p_adjL = multipletests(pL, method="fdr_bh")[1]
            except ZeroDivisionError:
                print >>sys.stderr, pL
                print >>sys.stderr, tmpL
                continue
            #---------------------------------------------
            for eachtmpL, p_adj in zip(tmpL, p_adjL):
                eachtmpL[-2] = format(eachtmpL[-2],'0.2E')
                eachtmpL[-1] = format(p_adj,'0.2E')
            #------------------------------------------------
            tmpL.sort(key=lambda x: float(x[-2]))
            
            headerL = [annoType + "_id", annoType + "_description"]
            for i in tmpL:
                if float(i[-1]) > pvalue_thresh:
                    break
                if headerL:
                    print >>fh_out, "%s\t%s\tTargetGene\tTargetCount\tTargetTotal\tTermCount\tTotalAnno\tOdds_ratio\tPvalue\tFDR" \
                        % (headerL[0], headerL[1])
                    headerL = []
                print >>fh_out, '\t'.join(i)

            #----------------------
            if out_file or grp != 1:
                fh_out.close()
                if topPlot:
                    os.system("head -n "+str(topPlot+1) + ' ' + out_file + '>' + out_file+'.top')
                    enrichPlot = ["sp_enrichmentPlot.sh -f ", out_file+'.top', '-o Odds_ratio -T numeric -v '+annoType+'_description -c FDR -s TargetCount -l FDR -x "Odds ratio" -y '+annoType]
                    os.system(' '.join(enrichPlot))
                    os.system(' '.join(['/bin/mv -f', out_file+'.top.scatterplot.dv.pdf', out_file[:-3]+'scatterplot.dv.pdf']))
            #--------------------------------------
        #-----------------------------------
    #--------------------------------------------------------
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


