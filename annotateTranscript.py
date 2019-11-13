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
    This is designed to annotate transcripts.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from tools import *

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
    parser.add_option("-o", "--output-prefix", dest="prefix",
        help="Prefix of the project. Any string contains numbers, alphabets and underlines only.")
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="Transcript FASTA file")
    parser.add_option("-p", "--pep-file", dest="pep",
        help="Translated fasta file. Optional.")
    parser.add_option("-s", "--strand-specific", dest="strand_specific",
        default=False, action="store_true",  help="Default both forward and reverse strand of given transcripts will be ananlyzed. If specified, only the forward strand will be analyzed.")
    parser.add_option("-g", "--gene_trans_map", dest="gene_trans_map",
        help="A two column file with first column as gene names and second column as transcript names. Optional. If not given, the sytem will treat each gene as one transcript.")
    parser.add_option("-m", "--mail", dest="mail",
        default="vip@ehbio.com", 
        help="DEfault <vip@ehbio.com>. Multiple mails can be separated by ','.")
    parser.add_option("-G", "--getGO", dest="GO",
        default=False, action="store_true", 
        help="Specify to get gene ontology information in separate file for GO enrichment analysis.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def outputMakefile(transcript, prefix, pep, strand_specific, gene_trans_map, go):
    if not pep:
        pep = "$(transcript).transdecoder.pep"
        os.system("touch FormatInputPep")
    else:
        os.system("touch Unigene_translation")
    if strand_specific:
        transdecoder_s = '-S'
    else:
        transdecoder_s = ''
    if gene_trans_map:
        os.system("touch Generate_gene_trans_map")
    else:
        gene_trans_map = "gene_trans_map"

    Makefile = "Makefile"
    Makefile_fh = open(Makefile, 'a')

    print >>Makefile_fh, '''
include /MPATHB/pipeline/makefile.tm

ifeq (makefile.am.template, $(wildcard makefile.am.template))
include makefile.am.template
endif

include /MPATHB/pipeline/Makefile.NGS.basic.NEW

prefix={prefix}
transcript={transcript}
final_fasta=$(transcript)
# strand-specific (only analyzes top strand)
# transdecoder_s=-S
transdecoder_s={transdecoder_s}
pep={pep}
final_pep=$(pep)
gene_trans_map?={gene_trans_map}

time:=$$(date +'%Y%m%d-%H%M%S')

Generate_gene_trans_map:
	touch $@
	grep '>' $(transcript) | sed 's/>//' | awk 'BEGIN{{OFS="\\t"}}{{print $$1, $$1}}' >$(gene_trans_map)

FormatInputPep:
	touch $@
	#pep file should contain `len:pep_length` and `scaffold:1-100(+)`
	# Only needed when pep file is not generated by TransDecoder
	/bin/mv $(pep) $(pep).$(time)
	formatFastaAddLen.py -i $(pep).$(time) -a >$(pep) 

Unigene_translation:
	touch $@
	TransDecoder -t $(transcript) -m 100 $(transdecoder_s) --search_pfam $(pfam) --CPU 20
	parseTransdecoder.pep2.py -f $(transcript).transdecoder.pep -m $(gene_trans_map) -s 1
	parseTransdecoder.pep2.py -f $(transcript).transdecoder.cds -m $(gene_trans_map)
	countTransdecoder.pep.py -i $(transcript).transdecoder.pep 
	$(mail)

Annotate_to_SwissProt_using_Blastx sprot_x:
	touch $@
	@echo "Begin $@ `date`"
	mkdir -p Trinotate
	blastx -query $(transcript) -task 'blastx-fast' -db $(uniprot_sprot) -num_threads 50 -max_target_seqs 1 -outfmt 6 -evalue 0.001 >Trinotate/$(prefix).trinity.uniprot_sprot.blastx 2>>Trinotate_$@.log
	@echo "End $@ `date`"
	$(mail)

Annotate_to_SwissProt_using_Blastp sprot_p:
	touch $@
	@echo "Begin $@ `date`"
	mkdir -p Trinotate
	blastp -query $(pep) -task 'blastp-fast' -db $(uniprot_sprot) -num_threads 20 -max_target_seqs 1 -outfmt 6 -evalue 0.001 >Trinotate/$(prefix).transdecoder.uniprot_sprot.blastp 2>>Trinotate_$@.log
	@echo "End $@ `date`"
	$(mail)


Annotate_to_TrEMBL_using_Blastx uniref90_x:
	touch $@
	@echo "Begin $@ `date`"
	mkdir -p Trinotate
	#blastx -query $(final_fasta) -task 'blastx-fast' -db $(uniprot_uniref90) -num_threads 30 -max_target_seqs 1 -outfmt 6 -evalue 0.00001 >Trinotate/$(prefix).trinity.uniprot_uniref90.blastx 2>>Trinotate_$@.log &
	touch Trinotate/$(prefix).trinity.uniprot_uniref90.blastx
	@echo "End $@ `date`"
	$(mail)

Annotate_to_TrEMBL_using_Blastp uniref90_p:
	touch $@
	@echo "Begin $@ `date`"
	mkdir -p Trinotate
	blastp -query $(final_pep) -task 'blastp-fast' -db $(uniprot_uniref90) -num_threads 20 -max_target_seqs 1 -outfmt 6 -evalue 0.0001 >Trinotate/$(prefix).transdecoder.uniprot_uniref90.blastp 
	@echo "End $@ `date`"
	$(mail)


Annotate_to_Pfam domain:
	touch $@
	@echo "Begin $@ `date`"
	mkdir -p Trinotate
	hmmscan --cpu 18 --domtblout Trinotate/$(prefix).transdecoder.pfam $(pfam) $(final_pep) >Trinotate_$@.log 2>&1
	@echo "End $@ `date`"
	$(mail)

signalp:
	touch $@
	@echo "Begin $@ `date`"
	mkdir -p Trinotate
	signalp -f short -n Trinotate/$(prefix).transdecoder.$@ $(final_pep) >Trinotate_$@.log 2>&1
	@echo "End $@ `date`"
	$(mail)

transmembrane:
	touch $@
	@echo "Begin $@ `date`"
	mkdir -p Trinotate
	tmhmm --short $(final_pep) >Trinotate/$(prefix).transdecoder.tmhmm 2>Trinotate_$@.log
	@echo "End $@ `date`"
	$(mail)


rRNA:
	touch $@
	@echo "Begin $@ `date`"
	mkdir -p Trinotate
	RnammerTranscriptome.pl --transcriptome $(final_fasta) --path_to_rnammer `which rnammer` >Trinotate_$@.log 2>&1
	/bin/mv *.rnammer.gff Trinotate/
	@echo "End $@ `date`"
	$(mail)

Trinotate_program=/MPATHB/soft/Trinotate-3.0.1/Trinotate
Trinotate_home=/MPATHB/soft/Trinotate-3.0.1/
load=yes
report=yes
go={go}

Trinotate Summarize_annotation:
	@echo "Begin $@ `date`"
	mkdir -p Trinotate
	touch $@
ifeq ($(load), yes)
	/bin/cp -f $(trinotate_sqlite) Trinotate/
	$(Trinotate_program) Trinotate/Trinotate.sqlite init --gene_trans_map $(gene_trans_map) --transcript_fasta $(transcript) --transdecoder_pep $(pep)
	if test -s Trinotate/$(transcript).rnammer.gff; then $(Trinotate_program) Trinotate/Trinotate.sqlite LOAD_rnammer Trinotate/$(transcript).rnammer.gff; fi
	$(Trinotate_program) Trinotate/Trinotate.sqlite LOAD_swissprot_blastp Trinotate/$(prefix).transdecoder.uniprot_sprot.blastp
	$(Trinotate_program) Trinotate/Trinotate.sqlite LOAD_swissprot_blastx Trinotate/$(prefix).trinity.uniprot_sprot.blastx
	$(Trinotate_program) Trinotate/Trinotate.sqlite LOAD_custom_blast --outfmt6 Trinotate/$(prefix).transdecoder.uniprot_uniref90.blastp --prog blastp --dbtype TrEMBL
	$(Trinotate_program) Trinotate/Trinotate.sqlite LOAD_custom_blast --outfmt6 Trinotate/$(prefix).trinity.uniprot_uniref90.blastx --prog blastx --dbtype TrEMBL
	#$(Trinotate_program) Trinotate/Trinotate.sqlite LOAD_trembl_blastp Trinotate/$(prefix).transdecoder.uniprot_uniref90.blastp
	#$(Trinotate_program) Trinotate/Trinotate.sqlite LOAD_trembl_blastx Trinotate/$(prefix).trinity.uniprot_uniref90.blastx
	$(Trinotate_program) Trinotate/Trinotate.sqlite LOAD_pfam Trinotate/$(prefix).transdecoder.pfam
	if test -s Trinotate/$(prefix).transdecoder.tmhmm; then $(Trinotate_program) Trinotate/Trinotate.sqlite LOAD_tmhmm Trinotate/$(prefix).transdecoder.tmhmm; fi
	if test -s Trinotate/$(prefix).transdecoder.signalp; then $(Trinotate_program) Trinotate/Trinotate.sqlite LOAD_signalp Trinotate/$(prefix).transdecoder.signalp; fi
endif
ifeq ($(report), yes)
	$(Trinotate_program) Trinotate/Trinotate.sqlite report >Trinotate/Trinotate_annotation_report.xls
	/bin/cp -f Trinotate/Trinotate_annotation_report.xls Trinotate/Trinotate_annotation_report.xls.original
	#--------------
	simplifyTrinotateAnnotation.py -i Trinotate/Trinotate_annotation_report.xls | awk 'BEGIN{{OFS=FS="\\t"}}ARGIND==1{{uni2g[$$1]=$$2;}}ARGIND==2{{if(FNR==1) $$2=$$2"\\tGene_name"; else {{name=uni2g[$$3]; if(name=="") {{split($$3,a,"_"); name=a[1];}} $$2=$$2"\\t"name;}} print $$0;}}' /MPATHB/resource/Trinotate/uniprot_sprot.id2gene - >Trinotate/Trinotate_annotation_report.simplify.xls
	s-plot vennDiagram -f Trinotate/Trinotate_annotation_report.xls.sta.xls -a TrEMBL -b Pfam -c KEGG -d SwissProt
	awk 'BEGIN{{OFS=FS="\\t"}}ARGIND==1{{iso2g[$$2]=$$1;}}ARGIND==2{{if($$2=="complete") id[iso2g[$$1]]=1;}}ARGIND==3{{if(id[$$1]==1) print $$0;}}' $(gene_trans_map) $(transcript).transdecoder.pep.id.xls Trinotate/Trinotate_annotation_report.xls.sta.xls >Trinotate/Trinotate_annotation_report.xls.complete.sta.xls
	s-plot vennDiagram -f Trinotate/Trinotate_annotation_report.xls.complete.sta.xls -a TrEMBL -b Pfam -c KEGG -d SwissProt
	tfFamily.py -i Trinotate/Trinotate_annotation_report.simplify.xls -t trinotate -f /MPATHB/resource/TFs/plant/PlantTFDB.family.std | sed 's/ /_/g' >Trinotate/Trinotate_annotation_report.TFfamily.xls
	awk 'BEGIN{{OFS=FS="\\t"}}{{if(FNR>1 && gene[$$1]=="") {{family[$$3]+=1; gene[$$1]=1;}}}}END{{print "TFfamily\\tCount"; for (tf in family) print tf,family[tf] | "sort -k2nr";}}' Trinotate/Trinotate_annotation_report.TFfamily.xls >Trinotate/Trinotate_annotation_report.TFfamily.sta.xls
	s-plot barPlot -f Trinotate/Trinotate_annotation_report.TFfamily.sta.xls -R 90 -x "TF families" -y "Number of TFs in each family" -P 'c(0.9, 0.9)'
endif
ifeq ($(go), True)
	$(Trinotate_home)/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls Trinotate/Trinotate_annotation_report.xls -G --include_ancestral_terms >Trinotate/Trinotate.gene.GO
	$(Trinotate_home)/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls Trinotate/Trinotate_annotation_report.xls -T --include_ancestral_terms >Trinotate/Trinotate.isoform.GO
	goPlotFromTrinotate.py -i Trinotate/Trinotate.gene.GO
endif
	@echo "END $@ `date`"
	$(mail)
'''.format(prefix=prefix, transcript=transcript, transdecoder_s=transdecoder_s, pep=pep, gene_trans_map=gene_trans_map, go=go)

    Makefile_fh.close()
#-----------outputMakefile----------------------------------

def outputAirflow(prefix, mail):
    airflow = prefix + '_trinotate_airflow.py'
    airflow_fh = open(airflow, 'w')
    print >>airflow_fh, '''
from airflow.models import DAG

from airflow.operators.bash_operator import BashOperator
from airflow.operators.email_operator import EmailOperator

from datetime import datetime, timedelta


one_min_ago = datetime.combine(datetime.today() - timedelta(minutes=1),
                                  datetime.min.time())
#now = datetime.now()

mail = {mail}
prefix = '{prefix}'
DAG_id = "Trinotate_" + prefix
dir="{dir}"

default_args = {{
    'owner': 'ct',         
    'depends_on_past': False, 
    'start_date': one_min_ago, 
    'email': [mail],
    'email_on_failure': True, 
    'email_on_retry': True, 
    'retries': 10, 
    'retry_delay': timedelta(hours=30), 
    #'queue': 'bash_queue',
    #'pool': 'backfill', 
    #'priority_weight': 10, 
    #'end_date': datetime(2016, 5, 29, 11, 30), 
}}

dag = DAG(DAG_id, default_args=default_args,
    schedule_interval="@once")


#Annotation---------------------------


Summarize_annotation = BashOperator(
    task_id="Summarize_annotation", 
    bash_command='(cd %s; make Summarize_annotation ) ' % dir, 
    dag=dag)

Generate_gene_trans_map = BashOperator(
    task_id="Generate_gene_trans_map", 
    bash_command='(cd %s; make Generate_gene_trans_map ) ' % dir, 
    dag=dag)
Generate_gene_trans_map.set_downstream([Summarize_annotation])


FormatInputPep = BashOperator(
    task_id="FormatInputPep", 
    bash_command='(cd %s; make FormatInputPep) ' % dir, 
    dag=dag)

Unigene_translation = BashOperator(
    task_id="Unigene_translation", 
    bash_command='(cd %s; make Unigene_translation) ' % dir, 
    dag=dag)

Annotate_to_SwissProt_using_Blastx =  BashOperator(
    task_id="Annotate_to_SwissProt_using_Blastx", 
    bash_command='(cd %s; make Annotate_to_SwissProt_using_Blastx) ' % dir, 
    dag=dag)

Annotate_to_SwissProt_using_Blastx.set_downstream([Summarize_annotation])

Annotate_to_SwissProt_using_Blastp =  BashOperator(
    task_id="Annotate_to_SwissProt_using_Blastp", 
    bash_command='(cd %s; make Annotate_to_SwissProt_using_Blastp) ' % dir, 
    dag=dag)

Annotate_to_SwissProt_using_Blastp.set_upstream([Unigene_translation, FormatInputPep])
Annotate_to_SwissProt_using_Blastp.set_downstream(Summarize_annotation)

Annotate_to_TrEMBL_using_Blastx =  BashOperator(
    task_id="Annotate_to_TrEMBL_using_Blastx", 
    bash_command='(cd %s; make Annotate_to_TrEMBL_using_Blastx) ' % dir, 
    dag=dag)

Annotate_to_TrEMBL_using_Blastx.set_downstream(Summarize_annotation)

Annotate_to_TrEMBL_using_Blastp =  BashOperator(
    task_id="Annotate_to_TrEMBL_using_Blastp", 
    bash_command='(cd %s; make Annotate_to_TrEMBL_using_Blastp) ' % dir, 
    dag=dag)

Annotate_to_TrEMBL_using_Blastp.set_upstream([FormatInputPep, Unigene_translation])
Annotate_to_TrEMBL_using_Blastp.set_downstream(Summarize_annotation)

Annotate_to_Pfam  = BashOperator(
    task_id="Annotate_to_Pfam", 
    bash_command='(cd %s; make Annotate_to_Pfam) ' % dir, 
    dag=dag)

Annotate_to_Pfam.set_upstream([FormatInputPep, Unigene_translation])
Annotate_to_Pfam.set_downstream(Summarize_annotation)

#---------------------------

Success_mail = EmailOperator(
    task_id="Success_mail",
    to=mail, 
    subject="%s Finished" % prefix, 
    html_content="{dir}", 
    dag=dag)

Success_mail.set_upstream(Summarize_annotation)

'''.format(dir=os.getcwd(), mail=mail, prefix=prefix)

    airflow_fh.close()
#------------outputAirflow--------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    transcript = options.filein
    prefix     = options.prefix
    pep        = options.pep
    strand_specific = options.strand_specific
    gene_trans_map  = options.gene_trans_map
    go         = options.GO
    mail       = options.mail.split(',')
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    checkLegalWord(prefix)

    outputMakefile(transcript, prefix, pep, strand_specific, gene_trans_map, go)
    outputAirflow(prefix, mail)
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


