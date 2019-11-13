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
    This is designed to get gene symbols for specific GO ids.

<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "rnorvegicus_gene_ensembl" interface = "default" >
		<Filter name = "go_parent_term" value = "GO:0015629,GO:0031532,GO:0005856"/>
		<Attribute name = "external_gene_name" />
		<Attribute name = "go_id" />
		<Attribute name = "name_1006" />
	</Dataset>
</Query>

http://www.ensembl.org/info/data/biomart/biomart_r_package.html

Dataset IDs:

"tguttata_gene_ensembl"	"Zebra Finch genes (taeGut3.2.4)"	"taeGut3.2.4"
"mmusculus_gene_ensembl"	"Mouse genes (GRCm38.p5)"	"GRCm38.p5"
"dmelanogaster_gene_ensembl"	"Fruitfly genes (BDGP6)"	"BDGP6"
"ptroglodytes_gene_ensembl"	"Chimpanzee genes (CHIMP2.1.4)"	"CHIMP2.1.4"
"hsapiens_gene_ensembl"	"Human genes (GRCh38.p10)"	"GRCh38.p10"
"celegans_gene_ensembl"	"Caenorhabditis elegans genes (WBcel235)"	"WBcel235"
"rnorvegicus_gene_ensembl"	"Rat genes (Rnor_6.0)"	"Rnor_6.0"


Filter IDs:

"name"	"description"
"1"	"chromosome_name"	"Chromosome/scaffold name"
"25"	"with_go"	"With GO ID(s)"
"61"	"ensembl_gene_id"	"Gene stable ID(s) [e.g. ENSG00000000003]"
"62"	"ensembl_transcript_id"	"Transcript stable ID(s) [e.g. ENST00000000233]"
"63"	"ensembl_peptide_id"	"Protein stable ID(s) [e.g. ENSP00000000233]"
"65"	"external_gene_name"	"Gene Name(s) [e.g. snoZ6]"
"66"	"external_transcript_name"	"Transcript Name(s) [e.g. snoZ6.3-201]"
"82"	"genedb"	"GeneDB ID(s) [e.g. Smp_009580.1:pep]"
"83"	"go"	"GO ID(s) [e.g. GO:0000002]"
"85"	"hgnc_symbol"	"HGNC symbol(s) [e.g. A1BG]"
"90"	"kegg_enzyme"	"KEGG Pathway and Enzyme ID(s) [e.g. 00010+1.1.1.1]"
"97"	"mirbase_id"	"miRBase ID(s) [e.g. hsa-let-7a-1]"
"100"	"entrezgene"	"NCBI gene ID(s) [e.g. 1]"
"199"	"go_parent_term"	"Parent term accession"
"200"	"go_parent_name"	"Parent term name"
"295"	"pfam"	"Pfam domain ID(s) [e.g. PF00001]"

Attribute IDs:

"ensembl_gene_id"	"Gene stable ID"	"feature_page"
"ensembl_transcript_id"	"Transcript stable ID"	"feature_page"
"description"	"Gene description"	"feature_page"
"chromosome_name"	"Chromosome/scaffold name"	"feature_page"
"start_position"	"Gene start (bp)"	"feature_page"
"end_position"	"Gene end (bp)"	"feature_page"
"strand"	"Strand"	"feature_page"
"transcript_start"	"Transcript start (bp)"	"feature_page"
"transcript_end"	"Transcript end (bp)"	"feature_page"
"transcription_start_site"	"Transcription start site (TSS)"	"feature_page"
"transcript_length"	"Transcript length (including UTRs and CDS)"	"feature_page"
"external_gene_name"	"Gene name"	"feature_page"
"percentage_gene_gc_content"	"% GC content"	"feature_page"
"go_id"	"GO term accession"	"feature_page"
"name_1006"	"GO term name"	"feature_page"
"definition_1006"	"GO term definition"	"feature_page"
"go"	"GO ID"	"feature_page"
"hgnc_symbol"	"HGNC symbol"	"feature_page"
"kegg_enzyme"	"KEGG Pathway and Enzyme ID"	"feature_page"
"entrezgene"	"NCBI gene ID"	"feature_page"
"pdb"	"PDB ID"	"feature_page"
"uniprotswissprot"	"UniProtKB/Swiss-Prot ID"	"feature_page"
"pfam"	"Pfam domain ID"	"feature_page"
"pfam_start"	"Pfam domain start"	"feature_page"
"pfam_end"	"Pfam domain end"	"feature_page"
"interpro"	"Interpro ID"	"feature_page"
"interpro_short_description"	"Interpro Short Description"	"feature_page"
"interpro_description"	"Interpro Description"	"feature_page"
"interpro_start"	"Interpro start"	"feature_page"
"interpro_end"	"Interpro end"	"feature_page"
"5_utr_start"	"5' UTR start"	"structure"
"5_utr_end"	"5' UTR end"	"structure"
"3_utr_start"	"3' UTR start"	"structure"
"3_utr_end"	"3' UTR end"	"structure"
"description"	"Gene description"	"structure"
"gene_biotype"	"Gene type"	"structure"
"exon_chrom_start"	"Exon region start (bp)"	"structure"
"exon_chrom_end"	"Exon region end (bp)"	"structure"
"5utr"	"5' UTR"	"sequences"
"3utr"	"3' UTR"	"sequences"
"gene_exon"	"Exon sequences"	"sequences"
"cdna"	"cDNA sequences"	"sequences"
"coding"	"Coding sequence"	"sequences"
"peptide"	"Peptide"	"sequences"
"upstream_flank"	"upstream_flank"	"sequences"
"downstream_flank"	"downstream_flank"	"sequences"
"ensembl_gene_id"	"Gene stable ID"	"sequences"
"description"	"Gene description"	"sequences"
"external_gene_name"	"Gene name"	"sequences"
"external_gene_source"	"Source of gene name"	"sequences"
"chromosome_name"	"Chromosome/scaffold name"	"sequences"
"start_position"	"Gene start (bp)"	"sequences"
"end_position"	"Gene end (bp)"	"sequences"
"gene_biotype"	"Gene type"	"sequences"
"family"	"Ensembl Protein Family ID(s)"	"sequences"
"uniparc"	"UniParc ID"	"sequences"
"uniprotswissprot"	"UniProtKB/Swiss-Prot ID"	"sequences"
"uniprotsptrembl"	"UniProtKB/TrEMBL ID"	"sequences"
"cdna_coding_start"	"CDS start (within cDNA)"	"sequences"
"cdna_coding_end"	"CDS end (within cDNA)"	"sequences"
"5_utr_start"	"5' UTR start"	"sequences"
"5_utr_end"	"5' UTR end"	"sequences"
"3_utr_start"	"3' UTR start"	"sequences"
"3_utr_end"	"3' UTR end"	"sequences"
"ensembl_transcript_id"	"Transcript stable ID"	"sequences"
"ensembl_peptide_id"	"Protein stable ID"	"sequences"
"transcript_biotype"	"Transcript type"	"sequences"
"strand"	"Strand"	"sequences"
"transcript_start"	"Transcript start (bp)"	"sequences"
"transcript_end"	"Transcript end (bp)"	"sequences"
"transcription_start_site"	"Transcription start site (TSS)"	"sequences"
"transcript_length"	"Transcript length (including UTRs and CDS)"	"sequences"
"cds_length"	"CDS Length"	"sequences"
"cds_start"	"CDS start"	"sequences"
"cds_end"	"CDS end"	"sequences"
"ensembl_exon_id"	"Exon stable ID"	"sequences"
"exon_chrom_start"	"Exon region start (bp)"	"sequences"
"exon_chrom_end"	"Exon region end (bp)"	"sequences"
"strand"	"Strand"	"sequences"
"rank"	"Exon rank in transcript"	"sequences"
"phase"	"Start phase"	"sequences"
"end_phase"	"End phase"	"sequences"
"cdna_coding_start"	"cDNA coding start"	"sequences"
"cdna_coding_end"	"cDNA coding end"	"sequences"
"genomic_coding_start"	"Genomic coding start"	"sequences"
"genomic_coding_end"	"Genomic coding end"	"sequences"
"is_constitutive"	"Constitutive exon"	"sequences"



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
    parser.add_option("-f", "--filter-name", dest="filter",
            metavar="FILEIN", help="Filter name and values separated by '@' like <go_parent_term@GO:0015629,GO:0031532,GO:0005856> (values are comma separated). Optional.")
    parser.add_option("-d", "--database", dest="database",
        help="Accept <rnorvegicus_gene_ensembl>, <mmusculus_gene_ensembl>, <hsapiens_gene_ensembl>")
    parser.add_option("-a", "--attributename", dest="attribute",
        help="Comma separated list of ENSEMBLE biomart supported names like <external_gene_name, go_id, name_1006>. See above for full lists and explanations.")
    parser.add_option("-o", "--output", dest="output",
        help="Output file name.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.database != None, "A database name needed for -d"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    filter = options.filter
    database = options.database
    attributeL = [i.strip() for i in options.attribute.split(',')]
    verbose = options.verbose
    output = options.output
    global debug
    debug = options.debug
    #-----------------------------------
    header = """http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "{}" interface = "default" >""".format(database)
    
    if filter:
        filterL = filter.split('@')
        filter_words = """<Filter name = "{}" value = "{}"/>""".format(filterL[0], filterL[1])
    else:
        filter_words = ""
    attribute_wordsL = ['<Attribute name = "{}" />'.format(attr) for attr in attributeL]
    footer = """</Dataset></Query>"""

    url = header+filter_words+''.join(attribute_wordsL)+footer
    #print >>sys.stderr, url
    
    cmd = 'wget -O '+output+' \''+url+'\''
    #print >>sys.stderr, cmd
    if os.system(cmd):
        print >>sys.stderr, cmd
        sys.exit(1)
        


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


