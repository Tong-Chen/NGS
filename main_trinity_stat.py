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
    This is designed to summarize the result of Trinity output.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
from tools import transferListToMultiLTable
#from multiprocessing.dummy import Pool as ThreadPool
from tools import *

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
    parser.add_option("-F", "--trinityFasta", dest="fasta",
        metavar="TRINITY-FASTA", help="The assembled transcripts. \
Normally <Unigene.fasta>. The program will generate \
'Unigene.fasta.transdecoder.pep' as Protein fasta file \
and 'Unigene.fasta.transdecoder.cds' as CDS fasta file.")
    parser.add_option("-i", "--stats", dest="filein",
        metavar="LENGTH", help="The output of TrinityStats.pl.")
    parser.add_option("-f", "--coverage", dest="fulllen",
        metavar="COMPLENETESS", help="The output of \
analyze_blastPlus_topHit_coverage.pl.")
    parser.add_option("-e", "--ExN50_result", dest="ExN50",
        help="Output of contig_E_statistic.pl")
    parser.add_option("-c", "--busco-result", dest="busco",
        help="Output of busco. Normally $(prefix).BUSCO.xls.")
#    parser.add_option("-g", "--gene_expr_estimate", dest="gene_expr",
#        help="Output of count_features_given_MIN_FPKM_threshold.pl \
#and sp_lines_lm.sh for genes. PNG file needed.")
#    parser.add_option("-t", "--transcript_expr_estimate", dest="iso_expr",
#        help="Output of count_features_given_MIN_FPKM_threshold.pl \
#and sp_lines_lm.sh for isoforms. PNG file needed.")
    parser.add_option("-S", "--sampleFile", dest="sampleFile",
        help="The sampleFile given to <PtR>.")
    parser.add_option("-s", "--sampleCor", dest="sampleCor",
        help="This parameter is used to get the heatmap and \
principle component analysis of given file. Normally \
<quantification/$(prefix).gene.expr.TMM.fpkm.matrix>.")
    parser.add_option("-a", "--anno", dest="anno",
        help="Single file or multiple files to be linked \
(separated by comma(,)). \
Normally \
<quantification/$(prefix).gene.expr.TMM.fpkm.matrix.anno.xls,\
quantification/$(prefix).isoform.expr.TMM.fpkm.matrix.anno.xls>.")
    parser.add_option("-A", "--append", dest="append",
        help="Description of files given to -a. \
Multiple descriptions are separated by `;`.")
    parser.add_option("-o", "--output-prefix", dest="out_prefix",
        help="The prefix of output files.")
    parser.add_option("-r", "--report-dir", dest="report_dir",
        default='report', help="Directory for report files. Default 'report'.")
    parser.add_option("-R", "--report-sub-dir", dest="report_sub_dir",
        default='2_assembling_quality', help="Directory for saving report figures and tables. This dir will put under <report_dir>,  so only dir name is needed. Default '2_assembling_quality'.")
    parser.add_option("-d", "--doc-only", dest="doc_only",
        default=False, action="store_true", help="Specify to only generate doc.")
    parser.add_option("-n", "--number", dest="number", type="int", 
        default=40, help="Set the maximum allowed samples for barplot. Default 40.\
 If more than this number of samples are given, heatmap will be used instead.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

#def copypng(dir, *fileL):
#    for pngFile in fileL:
#        if os.path.exists(pngFile):
#            os.system("cp -u %s %s" % (pngFile, dir))
#        if os.path.exists(pngFile.replace('png', 'pdf')):
#            os.system("cp -u %s %s" % (pngFile.replace('png', 'pdf'), dir))
##--------------------------------------------------------------
#
#def copypdf(dir, *fileL):
#    for pdfFile in fileL:
#        if os.path.exists(pdfFile):
#            os.system("cp -u %s %s" % (pdfFile, dir))
#            pngFile = pdfFile.replace('pdf', 'png')
#            if not os.path.exists(pngFile):
#                os.system(' '.join([
#                    "convert -density 150 -background white -alpha off",
#                    pdfFile,pngFile]))
#            os.system("cp -u %s %s" % (pngFile, dir))
#        else:
#            print >>sys.stderr, "Unexisted file %s" % pdfFile
##--------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    fasta = options.fasta
    file = options.filein
    ExN50 = options.ExN50
    #gene_expr = options.gene_expr
    #iso_expr = options.iso_expr
    busco = options.busco
    coverage = options.fulllen
    report_dir = options.report_dir
    sampleFileL = []
    for line in open(options.sampleFile):
        name = line.split('\t')[0]
        if name not in sampleFileL:
            sampleFileL.append(name)
    #------------------------------------------
    sampleCor = options.sampleCor
    annoL  = options.anno.split(',')
    appendL = options.append.split(';')
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    print '\n## 核心数据下载\n\n'

    dirname = "2_1_Unigene_transcript_assembl_fasta/"
    dir = report_dir+"/"+dirname
    os.system("mkdir -p %s" % dir)
    trinity = fasta.replace('Unigene', 'Trinity') 
    pep = fasta + ".transdecoder.pep"
    cds = fasta + ".transdecoder.cds"
    gs_pep = pep + ".gene_singleton.xls"
    is_pep = pep + ".isoform_singleton.xls"
    copyzip(dir, trinity, fasta, pep, cds, gs_pep, is_pep)
    
    print "为了方便下载，我们把组装的Unigene、蛋白、CDS的FASTA序列 (Table \@ref(tab:unigene-cds-prot-table)) 和样品基因表达注释表格 (Table \@ref(tab:expr-anno-total)) 放在了第一部分，后面再详细介绍拼装质量的评估、表达定量、差异分析等。"

    trinity2,fasta2,cds2,pep2,gs_pep2,is_pep2 = getRelativeDir([trinity,fasta,cds,pep,gs_pep,is_pep],dirname)

    print '''

Table: (\#tab:unigene-cds-prot-table) 组装的Unigene、蛋白、CDS的FASTA序列下载。

 ```{{r unigene-cds-prot-table, results="asis"}}
unigene_cds_prot_table <- "拼装基因序列下载;文件描述
[Unigene sequences]({fasta}.zip);FASTA[^1]格式的Unigene序列
[Unigene CDS sequences]({cds}.zip);FASTA[^2]格式的编码区序列
[Unigene protein sequences]({pep}.zip);FASTA格式蛋白序列
[Unigene protein sequences (excel)]({gs_pep}.zip);Excel[^3]表格形式的蛋白序列，方便使用Excel的vlookup函数批量提取序列
[All transcripts sequences]({trinity}.zip);FASTA[^1]格式的所有转录本序列"

unigene_cds_prot_table_mat <- read.table(text=unigene_cds_prot_table, sep=";", header=T)
knitr::kable(unigene_cds_prot_table_mat, booktabs=T, format="markdown")
```
'''.format(fasta=fasta2, cds=cds2, pep=pep2, gs_pep=gs_pep2, trinity=trinity2)

    #print "* [Unigene sequences (FASTA) FASTA格式的Unigene序列](%s%s.zip)\n" % \
    #    (dirname, os.path.split(fasta)[1])
    #print "* [Unigene CDS sequences (FASTA) FASTA格式的编码区序列](%s%s.zip)\n" % \
    #    (dirname, os.path.split(cds)[1])
    #print "* [Unigene protein sequences (FASTA) FASTA格式蛋白序列](%s%s.zip)\n" % \
    #    (dirname, os.path.split(pep)[1])
    #print "* [Unigene protein sequences (excel) Excel表格形式的蛋白序\
    #列，方便使用Excel的vlookup函数批量提取序列](%s%s.zip)\n" % \
    #    (dirname, os.path.split(gs_pep)[1])
    #print "* [Unigene protein sequences for all isoforms (excel)](%s%s.zip)\n" % \
    #    (dirname, os.path.split(is_pep)[1])
    
    print """
[^1]: FASTA格式的序列中，以 '>' 开头的行为序列的名字，其它行为序列，可直接用于Blast、ClustalW、Mega等分析软件。

[^2]: CDS预测的编码区序列，不包含UTR (untranslated regions)。

[^3]:如果把基因的表达量、注释和序列放在一个Excel表格中，文件会特别大，不容易打开。而如果序列用Fasta格式存储，又不方便提取；所以我们提供了Excel表格形式的序列，第一列为基因名字，第二列为对应序列，方便使用Excel的vlookup函数批量提取序列。当我们根据表达数据筛选到基因后就可以利用此文件批量获取基因的序列，保证了便利性。
"""

    dir_base = "2_5_Unigene_transcript_expression_annotation_profile/"
    
    dir = report_dir+"/%s" % dir_base
    os.system("mkdir -p %s" % dir)
    len_annoL = len(annoL)
    resultL_200 = ["File;Description"]
    for i in range(len_annoL):
        anno = annoL[i]
        append = appendL[i]
        copyzip(dir, anno)
        annoN = os.path.split(anno)[-1]
        tmp_anno = "[%s](%s.zip);%s" % (annoN, dir_base+annoN, append)
        resultL_200.append(tmp_anno)
    resultL_200 = '\n'.join(resultL_200)

    print """
    
Table: (\#tab:expr-anno-total) 所有样品的基因表达和注释的总表格。

```{r expr-anno-total, results="asis"}
expr_anno_total <- "%s"
expr_anno_total_mat <- read.table(text=expr_anno_total, sep=";", header=T)

knitr::kable(expr_anno_total_mat, booktabs=T, format="markdown")
```
""" % resultL_200

    print 
    print "\n## 重头组装序列质量评估\n"

    curation_label = "DE_novo_assemble_quality_estimation"
    knitr_read_txt(report_dir, curation_label)

    print "\n### 重头组装序列总体统计\n"
    
    print """
序列组装质量评估见 (Table \@ref(tab:assemble-quality-t))

Table: (\#tab:assemble-quality-t) Stats for all _de novo_ assembled transcript contigs. Given a set of contigs, each with its own length, the N50 length is defined as the length for which the collection of all contigs of that length or longer contains at least half of the sum of the lengths of all contigs, and for which the collection of all contigs of that length or shorter also contains at least 50% of the sum of the lengths of all contigs. N50常用于从序列长度的概念评估序列拼装的质量。其计算方式为把所有拼装出的序列按照从大到小排列，然后把其长度依次相加，直到获得的长度总和不低于所有序列长度总和一半时，所在转录本的长度定位N50。一般认为N50越长越好。但对于转录组拼装，获得的全长转录本的数目和转录本的注释率是评价拼装质量的更好的标准。

"""

    tableL = [["Item", "Value"]]
    for line in fh:
        if line.find("Stats based on") != -1:
            tableL.append([" * ", " * "])
        lineL = line.strip().split(':')
        if len(lineL) > 1:
            key, value = lineL
            key = key.strip()
            if key == "Percent GC":
                key = "GC content (%)"
            elif key.startswith("Contig"):
                key += " (nt)"
            elif key.startswith("Median"): 
                key += " (nt)"
            elif key.startswith("Average"):    
                key += " length (nt)"
            tableL.append([key, value.strip()])
        if line.find("Total assembled bases") != -1:
            break
    #-------------END reading file----------
    fill_star_1 = " "+"*" * 87+" "
    fill_star_2 = " "+"*" * 27+" "
    tableL.append([fill_star_1, fill_star_2])
    print '\n'.join(transferListToMultiLTable(tableL))
    print "     "
    
    dirname = "2_2_Unigene_transcript_assembl_quality_evaluation/"
    dir = report_dir+"/%s" % dirname

    os.system("mkdir -p %s" % dir)
    print 
    #if os.path.exists(en50):
    #    os.system("cp -u %s %s" % (en50, dir))
    #if os.path.exists(en50.replace('png', 'pdf')):
    #    os.system("cp -u %s %s" % (en50.replace('png', 'pdf'), dir))
    
    len_distrib = fasta + ".len_distrib.xls"
    len_distrib_name = os.path.split(len_distrib)[1]
    len_distrib_pdf = len_distrib + ".densityHist.pdf"
    len_distrib_pdf_name = os.path.split(len_distrib_pdf)[1]
    copypdf(dir, len_distrib_pdf)
    copy(dir, len_distrib)

    print "### 组装的UniGene的长度分布\n"

    print '''

(ref:unigene-length-distrib) Length distribution of Unigenes. The dashed line represents average length of all Unigenes. 直方图展示Unigene的长度分布；虚线表示所有Unigene的平均长度。[PDF](%s) [SOURCE](%s)

```{r unigene-length-distrib, fig.cap="(ref:unigene-length-distrib)"}
knitr::include_graphics("%s")
```
''' % (dirname+len_distrib_pdf_name, \
           dirname+len_distrib_name, \
           dirname+len_distrib_pdf_name.replace('pdf', 'png'))


    #en50 = ExN50 + "_N50.lines.pdf"
    #copypdf(dir, en50)
    #print "\n\n下面两幅图统计了不同表达水平的基因集的N50和绝对数目。左侧图结合基因的表达量，展示了表达最高的10%、20%、30%等基因的N50的数值，右侧的图展示的是每个对应的基因集内包含的基因数目。两者的结合关注的是对确实有可能发挥功能的基因的组装质量的评估。\n\n"
    #print "[ct-columns]\n"    
    #print "[ct-column=0.49]\n"    
    #print "[![N50 statistics for top most highly expressed transcripts.](%s)](%s)\n" \
    #    % (dirname + os.path.split(en50)[1].replace('pdf', 'png'), \
    #       dirname + os.path.split(en50)[1])
        
    #exprN50 = ExN50 + "_expr.lines.pdf"
    #copypdf(dir, exprN50)
    #if os.path.exists(expr):
    #    os.system("cp -u %s %s" % (expr, dir))
    #if os.path.exists(expr.replace('png', 'pdf')):
    #    os.system("cp -u %s %s" % (expr.replace('png', 'pdf'), dir))
    #print "[ct-column=0.49]\n"    
    #print "[![Cumulative counts of top most highly expressed \
#transcripts.](%s)](%s)\n" \
    #    % (dirname + os.path.split(exprN50)[1].replace('pdf', 'png'),\
    #        dirname + os.path.split(exprN50)[1])
    #print "[/ct-columns]\n"    
        
    print '### 组装的基因覆盖真核生物保守基因的比例\n'

    #busco_fh = open(busco)
    #line = busco_fh.readline()
    #while line.find("Statistics of the completeness of the genome") == -1:
    #    line = busco_fh.readline()
    #busco_fh.readline()
    #buscoL = []
    #headerL = ';'.join(busco_fh.readline().strip().split())
    #buscoL.append(headerL)
    #while 1:
    #    line = busco_fh.readline()
    #    line = line.strip()
    #    line = line.replace('Group ', 'Group_')
    #    if line.find("These results are based on") != -1:
    #        break
    #    if line:
    #        lineL = line.split()
    #        buscoL.append(';'.join(lineL))
    ##-------------------------------------------
    #busco_fh.close()
    #
    #busco_mat = '\n'.join(buscoL)

    busco_pdf = busco + ".stackBars.pdf"
    busco_pdf_name = os.path.split(busco_pdf)[1]

    copypdf(dir, busco_pdf)
    copy(dir, busco)

    #print >>sys.stderr, dirname, busco_pdf_name, busco

    print '''

(ref:busco-result) Quantitative assessment of genome assembly and annotation completeness based on evolutionarily informed expectations of gene content. The recovered genes are classified as ‘complete’ when their lengths are within two standard deviations of the BUSCO group mean length (i.e. within ∼95%% expectatio). ‘Complete’ genes found with more than one copy are classified as ‘duplicated’. These should be rare  as BUSCOs are evolving under single-copy control, and the recovery of many duplicates may therefore indicate erroneous assembly of haplotypes. Genes only partially recovered are classified as ‘fragmented’, and genes not recovered are classified as ‘missing’. Finally, the ‘number of genes used’ indicates the resolution and hence is informative of the confidence of these assessments. [PDF](%s) [SOURCE](%s)

```{r busco-result, fig.cap="(ref:busco-result)"}
knitr::include_graphics("%s")
```
''' % (dirname+busco_pdf_name, dirname + busco, 
           dirname+busco_pdf_name.replace('pdf', 'png'))

#    print '''
#```{{r busco-result}}
#busco_mat <- "{busco_mat}"
#
#busco_mat <- read.table(text=busco_mat, header=T, row.names=NULL, sep=';', check.names=F, comment="")
#knitr::kable(busco_mat, booktabs=T, caption="Statistics of the completeness of the assenbled Unigenes based on 248 CEGs. Complete represents Unigenes completely matched with CEGs. Partial represents Unigenes partially matched with CEGs. Prots = number of ultra-conserved CEGs present in assembled Unigenes; %Completeness = percentage of ultra-conserved CEGs present in assembled Unigenes relative to all 248 CEGs; Total = total number of CEGs present in assembled Unigenes including putative orthologs; Average = average number of orthologs per CEG ; %Ortho = percentage of detected CEGS that have more than 1 ortholog.")
#```
#
#'''.format(busco_mat=busco_mat)

    print '### 组装的基因的完整性和注释比例的评估\n'

    print """Table: (\#tab:unigene-completenes) Stats of the completeness of _de novo_ assembled Unigenes. 'Hit_pct_cov' indicates the coverage of known transcripts by _de novo_ assembled transcripts. 对组装的基因的完整性和注释比例的评估。'Hit_pct_cov'表示本项目拼装的基因与SwissProt数据库中已经注释的基因的匹配百分比。一般认为匹配度越高，越有可能拼装出的为全长序列。
"""

    table2 = [['Hit_pct_cov (%)', 'Count of transcripts', "Cumulative count of transcripts"]]
    for line in open(coverage):
        if line[0] == '#':
            continue
        percent, count, cumulative = line.split()
        percent = "%d-%s" % (int(percent)-10, percent)
        table2.append([percent, count, cumulative])
    table2.append(['*'*25, "*"*40, "*"*70])
    #--------END coverage------------------------
    print '\n'.join(transferListToMultiLTable(table2))
    print

    pep_sta = pep + ".sta.xls"
    pep_sta_name = os.path.split(pep_sta)[1]
    copy(dir, pep_sta)

#    print "[预测的UniGene编码的蛋白的完整性，`5prime_partial`表示蛋白N\
#端不完整, `internal`表示蛋白N端和C端都不完整，`3prime_partial`表示蛋白\
#C端不完整。`complete`表示完整的编码框。](%s)\n" % dir_name + pep_sta_name
    
    pep_sta_pdf = pep_sta + ".stackBars.pdf"
    pep_sta_pdf_name = os.path.split(pep_sta_pdf)[1]

    copypdf(dir, pep_sta_pdf)

    print '### 拼装基因的完整性评估\n'

    print '''

(ref:protein-completeness) Status of predicted proteins. 预测的UniGene编码的蛋白的完整性，`5prime_partial`表示蛋白N端不完整, `internal`表示蛋白N端和C端都不完整，`3prime_partial`表示蛋白C端不完整。`complete`表示完整的编码框。[PDF](%s) [SOURCE](%s)

```{r protein-completeness, fig.cap="(ref:protein-completeness)"}
knitr::include_graphics("%s")
```
''' % (dirname+pep_sta_pdf_name, dirname + pep_sta_name, 
           dirname+pep_sta_pdf_name.replace('pdf', 'png'))
    

    dir_base = "2_3_Unigene_transcript_annotation_profile/"
    
    dir = report_dir+"/%s" % dir_base

    os.system("mkdir -p %s" % dir)

    
    go_anno = "Trinotate/Trinotate.gene.GO.go_subparent.xls"
    go_pdf = go_anno+".barplot.pdf"
    venn_anno = "Trinotate/Trinotate_annotation_report.xls.sta.xls"
    venn_pdf = venn_anno+".vennDiagram.pdf"
    
    tf_anno = "Trinotate/Trinotate_annotation_report.TFfamily.sta.xls"
    tf_pdf = tf_anno + "barplot.pdf"
    copypdf(dir, go_pdf, tf_pdf, venn_pdf)
    copy(dir, go_anno, venn_anno, tf_anno)
    
    print "## 组装序列注释总结\n"


    curation_label = "DE_novo_assemble_annotation_estimation"
    knitr_read_txt(report_dir, curation_label)

    print "### Gene Ontology classification for all UniGenes\n"

    print """

(ref:go-sum) Summary of Gene Ontology annotation for all UniGenes. [PDF](%s) [SOURCE](%s)

```{r go-sum, fig.cap="(ref:go-sum)"}
knitr::include_graphics("%s")
```
""" % (dir_base+os.path.split(go_pdf)[1],
         dir_base+os.path.split(go_anno)[1],
         dir_base+os.path.split(go_pdf)[1].replace('pdf', 'png'))

    print "### Number of genes annotated to Swiss-Prot, TrEMBL, Pfam and KEGG\n"

    print """

(ref:swiss-trembl-pfam-sum) Summary of genes annotated to Swiss-Prot, TrEMBL and Pfam. [PDF](%s) [SOURCE](%s)

```{r swiss-trembl-pfam-sum, fig.cap="(ref:swiss-trembl-pfam-sum)"}
knitr::include_graphics("%s")
```
""" % (dir_base+os.path.split(venn_pdf)[1],
         dir_base+os.path.split(venn_anno)[1],
         dir_base+os.path.split(venn_pdf)[1].replace('pdf', 'png'))

    print "### Number of UniGenes in each TF family\n"

    print """

(ref:tf-family-sum) Summary of annotated transcription factor families. [PDF](%s) [SOURCE](%s)

```{r tf-family-sum, fig.cap="(ref:tf-family-sum)"}
knitr::include_graphics("%s")
```
""" % (dir_base+os.path.split(tf_pdf)[1],
         dir_base+os.path.split(tf_anno)[1],
         dir_base+os.path.split(tf_pdf)[1].replace('pdf', 'png'))


    dirname = "2_4_Unigene_quantification_and_sample_replication_estimation/"
    
    dir = report_dir+"/%s" % dirname

    os.system("mkdir -p %s" % dir)


#    copypng(dir, gene_expr, iso_expr)
#    print "## Counting numbers of expressed transcripts or genes\n"
#    
#    print "[ct-columns]\n"    
#    print "[ct-column=0.49]\n"    
#    print "![Number of expressed genes in given set.](%s)\n" % (dirname + os.path.split(gene_expr)[1])
#    print "[ct-column=0.49]\n"    
#    print "![Number of expressed transcritps in given set.](%s)\n" % (dirname + os.path.split(iso_expr)[1])
#    print "[/ct-columns]\n"    
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    
    print "\n## 样品重复性评估\n"

    curation_label = "DE_novo_assemble_sample_replication_estimation"
    knitr_read_txt(report_dir, curation_label)

    count = 0
    for name in sampleFileL:
        count += 1
        fourIn1_pdf = name + ".rep_compare.pdf"
        #os.system("touch "+" quantification/"fourIn1_png)
        exist = copypdf(dir, "quantification/"+fourIn1_pdf)
        if exist:

            bar_pdf = cor_pdf = ma_pdf = sca_pdf = fourIn1_pdf
            bar_png = name + ".rep_compare-0.png"
            cor_png = name + '.rep_compare-3.png'
            ma_png  = name + ".rep_compare-2.png"
            sca_png = name + ".rep_compare-1.png"

            copypng(dir, "quantification/"+bar_png, "quantification/"+cor_png, 
                    "quantification/"+ma_png, "quantification/"+sca_png)
        
            print "### %s replicates quality estimation\n" % name
            print '''

(ref:replicate-estimation-%d) 从左到右，从上到下，依次为: A. Statistics of mapped fragments. [PDF](%s); B. Pairwise MA plots. x-axis: mean log(counts-per-million). y-axis: log(fold_change). Data points more than 2-fold different are highlighted in red.[PDF](%s); C. Pearson correlation of replicates. [PDF](%s); D. Pairwise comparisons of replicates log(counts-per-million) values. Data points more than 2-fold different are highlighted in red. [PDF](%s)

```{r replicate-estimation-%d, out.width="49%s", fig.cap="(ref:replicate-estimation-%d)"}
knitr::include_graphics(c("%s", "%s", "%s", "%s"))
```
''' % (count, dirname + bar_pdf, dirname + ma_pdf, dirname + cor_pdf, 
        dirname + sca_pdf, count, '%', count, dirname + bar_png, dirname + ma_png, 
        dirname + cor_png, dirname + sca_png)
        #-------------------------------------------------------
    #-------END all samples-------------------------------------------------------
    quant_dat = sampleCor + ".bp.xls"
    quant_box = quant_dat + ".boxplot.noOutlier.pdf"
    copypdf(dir, quant_box)
    copy(dir, quant_dat)

    print "## 样品基因表达定量统计\n"
    print """

(ref:gene-expr-boxplot) Boxplot of gene expression in each sample. [PDF](%s) [SOURCE](%s)

```{r gene-expr-boxplot, fig.cap="(ref:gene-expr-boxplot)"}
knitr::include_graphics("%s")
```
""" % (dirname+os.path.split(quant_box)[1],
        dirname+os.path.split(quant_dat)[1],
        dirname+os.path.split(quant_box)[1].replace('pdf', 'png'))

    print "## 样品相似性评估和聚类分析\n"
    

    curation_label = "DE_novo_assemble_sample_cluster"
    knitr_read_txt(report_dir, curation_label)

    sampleCor_heat_dat = sampleCor + ".log2.sample_cor.dat"
    sampleCor_heat = sampleCor + ".log2.sample_cor_matrix.pdf"
    sampleCor_pc   = sampleCor + ".log2.prcomp.principal_components.pdf"
    copy(dir, sampleCor_heat_dat)
    copypdf(dir, sampleCor_heat)
    exist_pca = copypdf(dir, sampleCor_pc)
    

    print """

(ref:sample-correlation) Pearson correlation analysis of all samples. [PDF](%s) [SOURCE](%s)

```{r sample-correlation, fig.cap="(ref:sample-correlation)"}
knitr::include_graphics("%s")
```
""" % (dirname + os.path.split(sampleCor_heat)[1], 
        dirname + os.path.split(sampleCor_heat_dat)[1], 
        dirname + os.path.split(sampleCor_heat.replace("pdf","png"))[1]) 
            
    if exist_pca:
        print """

(ref:sample-correlation-pca) Cluster analysis of all samples using PCA (principle component analysis). [PDF](%s)

```{r sample-correlation-pca, fig.cap="(ref:sample-correlation-pca)"}
knitr::include_graphics("%s")
```
""" % (dirname + os.path.split(sampleCor_pc)[1], 
          dirname + os.path.split(sampleCor_pc.replace("pdf","png"))[1])

    
    
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


