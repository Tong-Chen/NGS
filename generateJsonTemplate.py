#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to generate pipeline configuration JSON file
    using given parameters.
'''

import sys
import os
from json import dumps as json_dumps
from json import loads as json_loads
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import re
#from multiprocessing.dummy import Pool as ThreadPool

#import readline is needed before importting rpy2.robjects
import readline 
#import rpy2.robjects as robjects
from textwrap import wrap as text_wrap

from itertools import combinations
from tools import *

#from bs4 import BeautifulSoup
reload(sys)
sys.setdefaultencoding('utf8')

debug = 0

makefileam = ""

ignoreKeyD = {'Attention':1}

# Sample specific parameters 
customizeD = {'seq_type': 1, 'libtype':1,
    'execute_cufflinks':1, "execute_htseq": 1, "peak":1} 



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
        metavar="FILEIN", help="JSON file")
    parser.add_option("-t", "--type", dest="hts_type",
        type="choice", 
        choices=["RNA_seq:reference_based", "RNA_seq:denovo", 
            "Exome_seq", "ChIP_seq", "ctDNA", "RNA_seq:scRNA", 
            "RRBS", "WBS"], 
        help="Specify the type of data for analysis. \
Currently ['RNA_seq:reference_based', 'RNA_seq:denovo', \
'Exome_seq', 'ChIP_seq', 'ctDNA', 'RNA_seq:scRNA', 'RRBS', 'WBS'] are supprted. The colon symbol \
is used to define subtypes.")
    parser.add_option("-A", "--assembl-only", dest="assembl_only",
        default=False, action="store_true", 
        help="Specify to run transcriptome assemble only.")
    parser.add_option("-p", "--process", dest="process",
        help="A JSON file to specify the jobs to be run \
and their runway.")
    parser.add_option("-c", "--no-check", dest="check",
        default=True, action="store_false", 
        help="Check if makefile.am.template exists. Default <check>. Specify this parameter to overide previous makefile.am.template")
    parser.add_option("-D", "--dryRun", dest="dryRun",
        default=False, action="store_true", help="Dry run")
    parser.add_option("-o", "--owner", dest="owner",
        default='ct', help="Supply user name")
    parser.add_option("-q", "--pipeline-path", dest="ppPath",
        default="/home/ct/pipeline/", help="Default /home/ct/pipeline/")
    parser.add_option("-a", "--alternative-splicing", dest="AS",
        default=False, action="store_true", help="Do alternative splicing analysis or not. \
Default False. Specify to open this operation. Please pay attention that \
you have sent the right <split_len> and understand that you may lose many reads \
information.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def getVersion(key, value):
    """ 
    Get software version by running related commands.

    Args:
        key(str): software name
        value(str): the way to get the version.
                * `bash commands`: A bash command which return the version of software
                * `R package`: Indicate this is an R package
                * `version 3.0.9`: Directly supply the version
            
    Returns:
        version(str)
            
    Raises:
        None
    """
    if value == "R package":
        version = 'Rscript -e ' + '\'installed.packages()[c("'+key+'"),c("Package","Version")]\''.format(key)
        #print >>sys.stderr, version
        versionL = os.popen(version)
        if versionL:
            version = list(versionL)[1].strip().replace('"', '')
        else:
            version = "Uninstall"
        #version = robjects.r('installed.packages()[c("%s"),\
        #    c("Package","Version")]' % key)[1]
    elif value.startswith('version '):
        version = value
    else:
        tmp143 = os.popen(value).readlines()
        if tmp143:
            version = tmp143[0].strip()
        else:
            version = 'Uninstall'
    return version
#------------------------------------------------------------
def parseLibType(saveD, libtype):
    """ 
    Generate strand information for all programs.

    Args:
        saveD(dict): key-pair representing the program and their strand parameter
        libtype(str): unstrand, RF, R, FR, F

    Returns:
        None
    """
    saveD['pasa_strand'] = libtype 
    if libtype == 'unstrand':
        saveD['tophat_libtype'] = 'fr-unstranded'
        saveD['trinity_strand'] = ''
        saveD['trinity_quantification_par'] = ''
        saveD['igv_cnt_strand'] = ''
        saveD['transDecoder_strand'] = ''
        saveD['htseq_strand']   = 'no' 
        saveD['outWigStrand']   = "Unstranded"
        saveD["outSAMstrandField"] = "intronMotif"
        saveD['tpm_strand'] = 'none'
    elif libtype == "RF" or libtype == 'R':
        saveD['tophat_libtype'] = 'fr-firststrand'
        saveD['trinity_strand'] = '--SS_lib_type '+libtype
        saveD['trinity_quantification_par'] = '--SS_lib_type '+libtype
        saveD['igv_cnt_strand'] = '--strands first'
        saveD['htseq_strand']   = 'reverse'
        saveD['correct'] = 'yes'
        saveD['transDecoder_strand'] = '-S'
        saveD['outWigStrand']   = "Stranded"
        saveD["outSAMstrandField"] = "None"
        saveD['tpm_strand'] = 'reverse'
    elif libtype == "FR" or libtype == 'F':
        saveD['tophat_libtype'] = 'fr-secondstrand'
        saveD['trinity_strand'] = '--SS_lib_type '+libtype
        saveD['trinity_quantification_par'] = '--SS_lib_type '+libtype
        saveD['igv_cnt_strand'] = '--strands first'
        saveD['htseq_strand']   = 'yes'
        saveD['correct'] = 'yes'
        saveD['transDecoder_strand'] = '-S'
        saveD['outWigStrand']   = "Stranded"
        saveD["outSAMstrandField"] = "None"
        saveD['tpm_strand'] = 'forward'
    #--------------------------------------
#----------------------------------------------

def dealWithSpeciesInformation(value):
    value = value.upper()
    assembleD = {'GRCH38':{"species_name":"hsa", "enrich_db":"org.Hs.eg.db", "macs_gsize":"2.7e9"}, 
            'MM10': {"species_name":"mmu", "enrich_db":"org.Mm.eg.db", "macs_gsize":"1.87e9"},    
            'TAIR10': {"species_name":"ath", "enrich_db":"org.At.tair.db", "macs_gsize":"1e8"},    
            'DM6': {"species_name":"dme", "enrich_db":"org.Dm.eg.db", "macs_gsize":"1.2e8"},    
            'CE11': {"species_name":"cel", "enrich_db":"org.Ce.eg.db", "macs_gsize":"9e7"},     
            'RNOR6': {"species_name":"rno", "enrich_db":"org.Rn.eg.db" }    
    }
    subD = assembleD.get(value)
    if subD:
        #print >>sys.stderr, subD
        for key, value in subD.items():
            print >>makefileam, "{}={}".format(key, value)
    else:
        print >>sys.stderr, "*** No species specific information recorded. ***"

#-------------------------------------------
def output_innerD(sec_key, value):
    print >>makefileam, "#", sec_key
    for key_tmp, value_tmp in value.items():
        if isinstance(value_tmp, list):
            value_tmp = ' '.join(value_tmp)
        print >>makefileam, "{}={}".format(key_tmp, value_tmp)

#---output_innerD----------------------------

def traverseDict(key, valueD, level, sampleCustomizeD, hts_typeL):
 
    """ 
    Trave a dict and print all key-value pair.

    Args:
        key(str): keyname
        valueD(dict): a Dict
        level(int): Indicate the descendants
        sampleCustomizeD(dict): {samp1:{libtype:'R'};
            samp2:{libtype:'F'}}
        hts_typeL(list):[Exome_seq] or [ChIP_seq] or [RNA_seq, reference_based]

    Returns:
        None
    """
    print >>makefileam, "\n%s %s\n" % ('#'*level, key)
    commentD = valueD.get('comment', {})
    for sec_key, value in valueD.items():
        #-------Select the analysis type--------
        if key == "Analysis_type" and sec_key not in hts_typeL:
            continue
        #-------Select analysis sub-type--------
        
        if sec_key in ["comment", "SRA_list", "Sample_list", "Compare_pair"]:
            continue
        if isinstance(value, dict):
            if key in hts_typeL: 
                if sec_key in ['Treat_control']:
                    output_innerD(sec_key, value)
                    continue
                if sec_key not in hts_typeL: 
                    continue
                else: 
                    hts_typeL.remove(sec_key)
            traverseDict(sec_key, value, level+1, sampleCustomizeD,
                hts_typeL)
            #print >>makefileam, "\n"
        else:
            if value == 'default':
                continue
            if sec_key == 'version':
                value = getVersion(key, value)
                sec_key = key + '_version'
            if isinstance(value, list): 
                #join multiple parameters for programs
                #print >>sys.stderr, value
                value = ' '.join(value)
            if sec_key in customizeD:
                for valueD in sampleCustomizeD.values():
                    if sec_key == 'libtype':
                        parseLibType(valueD, value)
                    valueD[sec_key] = value
            else:
                comment_v = text_wrap(commentD.get(sec_key,""),40)
                #comment_v = unicode(comment_v)
                print >>makefileam, "\n#--", "\n#-- ".join(comment_v)
                print >>makefileam, "%s=%s" % (sec_key, str(value))
                #--------Fill in assemble related information-------
                if sec_key == "assembl":
                    dealWithSpeciesInformation(value)
        #-----------------------------------------------------                
    #-----------------------------------------------------                
#-------------------------------------

def parseSampleInfo(sampleD, work_dir):
    """ 
    Parse sample information and generate one file:
        compare_pair
        #sampleFile, compare_pair, sampleFile_denovo

    Args:
        sampleD(dict): A specilized dict
            
    Returns:
        sampleS(list): A list of samples
            
    Raises:
        None

    sampleFile
    Samp    conditions
    T0_1    T0
    T0_2    T0
    T0_3    T0
    T2_1    T2
    T2_2    T2
    T2_3    T2

    compare_pair
    T0      T2
    T0      T4
    T0      T12
    T0      T24
    T0      T48

    sampleFile_denovo
    T0      T0_1
    T0      T0_2
    T0      T0_3
    T2      T2_1
    T2      T2_2
    T2      T2_3

    """ 
    sampleL = []
    replicates = 'no'
    sample_list = sampleD["Sample_list"]
    commentD = sampleD.get('comment', {})
    #sampleFile = open("sampleFile", 'w')
    compare_pair = open(work_dir+"/compare_pair", 'w')
    #sampleFile_denovo = open("sampleFile_denovo", "w")
    
    #print >>sampleFile, "Samp\tconditions"
    
    for innerD in sample_list:
        key, valueL = innerD.items()[0]
        sampleL.extend(valueL)
        #for value in valueL:
        #    print >>sampleFile, "%s\t%s" % (value, key)
        #    print >>sampleFile_denovo, "%s\t%s" % (key, value)
    #sampleFile.close()
    #sampleFile_denovo.close()
    
    print >>makefileam, "\n# Sample_info"
    print >>makefileam, "\n#--", \
        "\n#-- ".join(text_wrap(commentD.get('Sample_list', ''), 50))
    print >>makefileam, "nameL=%s" % ' '.join(sampleL) 

    if len(sampleL) > len(sample_list):
        print >>makefileam, "replicates=yes"
    else:
        print >>makefileam, "replicates=no"

    if len(sampleL) > 4:
        print >>makefileam, "prin_comp=yes"
    else:
        print >>makefileam, "prin_comp=no"

    de_analysis = 0
    compare_pairLL = sampleD.get('Compare_pair',[])
    for itemL in compare_pairLL:
        if len(itemL)>1:
            de_analysis = 1
            item_comb = list(combinations(itemL, 2))
            print >>compare_pair, "\n".join(["\t".join(item) for item in item_comb])
    compare_pair.close()
    if de_analysis:
        print >>makefileam, "\n# Perfom DE analysis as indicated in `compare_pair`"
        print >>makefileam, "DE_analysis=%d" % de_analysis
    return sampleL, sample_list
#-----------------------------------------

def sampleFileGeneration(sample_list, work_dir, sampleL=[]):
    """ 
    Parse sample information and generate two files:
        sampleFile, sampleFile_denovo

    Args:
        sampleD(dict): A specilized dict
            
    Returns:
        sampleS(list): A list of samples
            
    Raises:
        None

    sampleFile
    Samp    conditions
    T0_1    T0
    T0_2    T0
    T0_3    T0
    T2_1    T2
    T2_2    T2
    T2_3    T2

    sampleFile_denovo
    T0      T0_1
    T0      T0_2
    T0      T0_3
    T2      T2_1
    T2      T2_2
    T2      T2_3

    """ 
    sampleFile = open(work_dir+"/sampleFile", 'w')
    sampleFile_denovo = open(work_dir+"/sampleFile_denovo", "w")
    
    print >>sampleFile, "Samp\tconditions"

    for innerD in sample_list:
        key, valueL = innerD.items()[0]
        for value in valueL:
            if (sampleL and value in sampleL) or (not sampleL):
                assert key==value or value.find(key+'_') != -1, 'Wrong type '+value+'. Replicates must in format sampleName_rep.'
                print >>sampleFile, "%s\t%s" % (value, key)
                print >>sampleFile_denovo, "%s\t%s" % (key, value)
    sampleFile.close()
    sampleFile_denovo.close()
#==============================================================

def airflow_etl_rnaseq_denovo(prefix, customerMailLs, ehbioMailLs, work_dir, sampleLs, dryRun, owner, assembl_only):
    count = 1
    output = work_dir + '/' + prefix + "_" + str(count) + "_"+ owner +'_denovo_rnaseq_airflow_dag.py'
    if dryRun:
        makeadd = "-n"
    else:
        makeadd = ""
    while os.path.exists(output):
        print >>sys.stderr, output+" exists! Please check first."
        count += 1
        output = work_dir + '/' + prefix + "_" + str(count) + "_"+ owner \
                +'_denovo_rnaseq_airflow_dag.py'
        #sys.exit(1)

    output_fh = open(output,  'w')
    if assembl_only:
        print >>output_fh, '''   
from airflow.models import DAG

from airflow.operators.bash_operator import BashOperator
from airflow.operators.email_operator import EmailOperator

from datetime import datetime, timedelta


one_min_ago = datetime.combine(datetime.today() - timedelta(minutes=1),
                                  datetime.min.time())
#now = datetime.now()

mail = [{mail}]
cc_mail = [{customerMail}]
prefix = '{prefix}'
DAG_id = "De_novo_RNAseq_{count}_" + prefix + "_{owner}"

## Rewrite following tow lines
dir = "{work_dir}"
nameL = [{sampleL}]


default_args = {{
    'owner': '{owner}',         
    'depends_on_past': False, 
    'start_date': one_min_ago, 
    'email': mail,
    'email_on_failure': True, 
    'email_on_retry': True, 
    'retries': 1000, 
    'retry_delay': timedelta(minutes=300), 
    #'queue': 'bash_queue',
    #'pool': 'backfill', 
    #'priority_weight': 10, 
    #'end_date': datetime(2016, 5, 29, 11, 30), 
}}

dag = DAG(DAG_id, default_args=default_args,
    schedule_interval="@once")

trinity_cat_fastq = BashOperator(
    task_id="De_novo_assemble", 
    bash_command='(cd %s; make {makeadd} trinity.cat_fastq) ' % dir, 
    retry_delay=timedelta(minutes=180),
    retries=1000, 
    dag=dag)


for i in nameL:
    SRA_download = BashOperator(
        task_id="Download_SRA_data_"+i, 
        bash_command='(cd %s; make {makeadd} %s.SRA_download sample=%s) ' % (dir, i, i), 
        retry_delay=timedelta(minutes=1),
        retries=100, 
        dag=dag)

    SRA_rename = BashOperator(
        task_id="Rename_SRA_data_"+i, 
        bash_command='(cd %s; make {makeadd} %s.SRA_rename sample=%s) ' % (dir, i, i), 
        dag=dag)

    SRA_rename.set_upstream(SRA_download)

    fastqc = BashOperator(
        task_id="Quality_estimation_for_sample_"+i, 
        bash_command='(cd %s; make {makeadd} %s.fastqc sample=%s) ' \
            % (dir, i, i), 
        dag=dag)

    fastqc.set_upstream(SRA_rename)

    preprocess = BashOperator(
        task_id="Trim_adaptors_lowquality_bases_for_sample_"+i, 
        bash_command='(cd %s; make {makeadd} %s.trim_noref sample=%s) ' \
            % (dir, i, i), 
        dag=dag)

    preprocess.set_upstream(fastqc)
    preprocess.set_downstream([trinity_cat_fastq])


#------------quantification---------------------------

Success_mail = EmailOperator(
    task_id="Success_mail",
    to=cc_mail, 
    cc=mail, 
    subject="%s Finished" % prefix, 
    html_content='Dear Sir/Mandam<br> Please check your results in <a href="http://210.74.4.67:11521/result/{prefix}_EHBio">http://210.74.4.67:11521/result/{prefix}_EHBio</a> (username: ehbio; password: {prefix}_EHBio) <br> EHBIO Gene Technology'.format(prefix=prefix), 
    dag=dag)

Success_mail.set_upstream(trinity_cat_fastq)
'''.format(mail=ehbioMailLs, customerMail=customerMailLs, prefix=prefix, 
        work_dir=work_dir, sampleL=sampleLs, makeadd=makeadd, owner=owner, 
        count = count)

    else:
        print >>output_fh, '''   
from airflow.models import DAG

from airflow.operators.bash_operator import BashOperator
from airflow.operators.email_operator import EmailOperator

from datetime import datetime, timedelta


one_min_ago = datetime.combine(datetime.today() - timedelta(minutes=1),
                                  datetime.min.time())
#now = datetime.now()

mail = [{mail}]
cc_mail = [{customerMail}]
prefix = '{prefix}'
DAG_id = "De_novo_RNAseq_{count}_" + prefix + "_{owner}"

## Rewrite following tow lines
dir = "{work_dir}"
nameL = [{sampleL}]


default_args = {{
    'owner': '{owner}',         
    'depends_on_past': False, 
    'start_date': one_min_ago, 
    'email': mail,
    'email_on_failure': True, 
    'email_on_retry': True, 
    'retries': 1000, 
    'retry_delay': timedelta(minutes=300), 
    #'queue': 'bash_queue',
    #'pool': 'backfill', 
    #'priority_weight': 10, 
    #'end_date': datetime(2016, 5, 29, 11, 30), 
}}

dag = DAG(DAG_id, default_args=default_args,
    schedule_interval="@once")

Doc_all = BashOperator(
    task_id="Doc_all",
    bash_command='(cd %s; make {makeadd} DOC) ' % dir,
    dag=dag) 

Doc_sequencing_quality = BashOperator(
    task_id="Doc_sequencing_quality", 
    bash_command='(cd %s; make {makeadd} a_main_fastqc_doc) ' % dir, 
    dag=dag)

Doc_sequencing_quality.set_downstream(Doc_all)

Doc_assemble_quality = BashOperator(
    task_id="Doc_assemble_quality", 
    bash_command='(cd %s; make {makeadd} b_denovo_assembl_quality b_denovo_assembl_quality_doc) ' % dir, 
    dag=dag)

Doc_assemble_quality.set_downstream(Doc_all)


trinity_cat_fastq = BashOperator(
    task_id="De_novo_assemble", 
    bash_command='(cd %s; make {makeadd} trinity.cat_fastq) ' % dir, 
    retry_delay=timedelta(minutes=180),
    retries=1000, 
    dag=dag)


for i in nameL:
    SRA_download = BashOperator(
        task_id="Download_SRA_data_"+i, 
        bash_command='(cd %s; make {makeadd} %s.SRA_download sample=%s) ' % (dir, i, i), 
        retry_delay=timedelta(minutes=1),
        retries=100, 
        dag=dag)

    SRA_rename = BashOperator(
        task_id="Rename_SRA_data_"+i, 
        bash_command='(cd %s; make {makeadd} %s.SRA_rename sample=%s) ' % (dir, i, i), 
        dag=dag)

    SRA_rename.set_upstream(SRA_download)

    Statistics_sequenced_reads_bases = BashOperator(
        task_id="Statistics_sequenced_reads_bases_"+i, 
        bash_command='(cd %s; make {makeadd} %s.statistics_sequenced_reads_bases sample=%s) ' % (dir, i, i), 
        dag=dag)

    Statistics_sequenced_reads_bases.set_upstream(SRA_rename)
    Statistics_sequenced_reads_bases.set_downstream(Doc_sequencing_quality)

    fastqc = BashOperator(
        task_id="Quality_estimation_for_sample_"+i, 
        bash_command='(cd %s; make {makeadd} %s.fastqc sample=%s) ' \
            % (dir, i, i), 
        dag=dag)

    fastqc.set_upstream(SRA_rename)

    preprocess = BashOperator(
        task_id="Trim_adaptors_lowquality_bases_for_sample_"+i, 
        bash_command='(cd %s; make {makeadd} %s.trim_noref sample=%s) ' \
            % (dir, i, i), 
        dag=dag)

    preprocess.set_upstream(fastqc)
    preprocess.set_downstream([trinity_cat_fastq, Doc_sequencing_quality])


#------------quantification---------------------------
Prepare_reference_transcriptome = BashOperator(
    task_id="Prepare_reference_transcriptome", 
    bash_command='(cd %s; make {makeadd} Prepare_reference_transcriptome) ' % dir, 
    dag=dag)

Prepare_reference_transcriptome.set_upstream(trinity_cat_fastq)

Integrate_multiple_samples_expression = BashOperator(
    task_id="Integrate_multiple_samples_expression", 
    bash_command='(cd %s; make {makeadd} Integrate_multiple_samples_expression) ' % dir, 
    dag=dag)

for i in nameL:
    task1 = BashOperator(
        task_id="Quantification_for_sample_"+i, 
        bash_command='(cd %s; make {makeadd} %s.trinity_quantification sample=%s) ' \
            % (dir, i, i), 
        dag=dag)
    task1.set_upstream(Prepare_reference_transcriptome)
    task1.set_downstream(Integrate_multiple_samples_expression)


Filter_lowly_expressed_Unigenes = BashOperator(
    task_id="Filter_lowly_expressed_Unigenes", 
    bash_command='(cd %s; make {makeadd} Filter_lowly_expressed_Unigenes) ' % dir, 
    dag=dag)

Filter_lowly_expressed_Unigenes.set_upstream(Integrate_multiple_samples_expression)

#Annotation---------------------------

Summarize_annotation = BashOperator(
    task_id="Summarize_annotation", 
    bash_command='(cd %s; make {makeadd} Summarize_annotation ) ' % dir, 
    dag=dag)

Unigene_translation = BashOperator(
    task_id="Unigene_translation", 
    bash_command='(cd %s; make {makeadd} Unigene_translation) ' % dir, 
    dag=dag)

Unigene_translation.set_upstream(Filter_lowly_expressed_Unigenes)

Summarize_SwissProt_coverage = BashOperator(
    task_id="Summarize_SwissProt_coverage", 
    bash_command='(cd %s; make {makeadd} Summarize_SwissProt_coverage) ' % dir, 
    dag=dag)

Summarize_SwissProt_coverage.set_downstream(Doc_assemble_quality)

Annotate_to_SwissProt_using_Blastx =  BashOperator(
    task_id="Annotate_to_SwissProt_using_Blastx", 
    bash_command='(cd %s; make {makeadd} Annotate_to_SwissProt_using_Blastx) ' % dir, 
    dag=dag)

Annotate_to_SwissProt_using_Blastx.set_upstream(Filter_lowly_expressed_Unigenes)
Annotate_to_SwissProt_using_Blastx.set_downstream([Summarize_annotation,
    Summarize_SwissProt_coverage])

Annotate_to_SwissProt_using_Blastp =  BashOperator(
    task_id="Annotate_to_SwissProt_using_Blastp", 
    bash_command='(cd %s; make {makeadd} Annotate_to_SwissProt_using_Blastp) ' % dir, 
    dag=dag)

Annotate_to_SwissProt_using_Blastp.set_upstream(Unigene_translation)
Annotate_to_SwissProt_using_Blastp.set_downstream(Summarize_annotation)

Annotate_to_TrEMBL_using_Blastx =  BashOperator(
    task_id="Annotate_to_TrEMBL_using_Blastx", 
    bash_command='(cd %s; make {makeadd} Annotate_to_TrEMBL_using_Blastx) ' % dir, 
    dag=dag)

Annotate_to_TrEMBL_using_Blastx.set_upstream(Filter_lowly_expressed_Unigenes)
Annotate_to_TrEMBL_using_Blastx.set_downstream(Summarize_annotation)

Annotate_to_TrEMBL_using_Blastp =  BashOperator(
    task_id="Annotate_to_TrEMBL_using_Blastp", 
    bash_command='(cd %s; make {makeadd} Annotate_to_TrEMBL_using_Blastp) ' % dir, 
    dag=dag)

Annotate_to_TrEMBL_using_Blastp.set_upstream(Unigene_translation)
Annotate_to_TrEMBL_using_Blastp.set_downstream(Summarize_annotation)

Annotate_to_Pfam  = BashOperator(
    task_id="Annotate_to_Pfam", 
    bash_command='(cd %s; make {makeadd} Annotate_to_Pfam) ' % dir, 
    dag=dag)

Annotate_to_Pfam.set_upstream(Unigene_translation)
Annotate_to_Pfam.set_downstream(Summarize_annotation)

#---------------------------
Evaluate_de_novo_assemble_quality = BashOperator(
    task_id="Evaluate_de_novo_assemble_quality", 
    bash_command='(cd %s; make {makeadd} Evaluate_de_novo_assemble_quality) ' % dir, 
    dag=dag)

Evaluate_de_novo_assemble_quality.set_upstream(Filter_lowly_expressed_Unigenes)
Evaluate_de_novo_assemble_quality.set_downstream(Doc_assemble_quality)

Core_Eukaryotic_Genes_Mapping_Approach = BashOperator(
    task_id="Core_Eukaryotic_Genes_Mapping_Approach", 
    bash_command='(cd %s; make {makeadd} Core_Eukaryotic_Genes_Mapping_Approach) ' % dir, 
    dag=dag)

Core_Eukaryotic_Genes_Mapping_Approach.set_upstream(Filter_lowly_expressed_Unigenes)
Core_Eukaryotic_Genes_Mapping_Approach.set_downstream(Doc_assemble_quality)


Evaluate_quantification_result = BashOperator(
    task_id="Evaluate_quantification_result", 
    bash_command='(cd %s; make {makeadd} Evaluate_quantification_result) ' % dir, 
    dag=dag)

Evaluate_quantification_result.set_upstream(Filter_lowly_expressed_Unigenes)
Evaluate_quantification_result.set_downstream(Doc_assemble_quality)

Evaluate_biological_reproducibility = BashOperator(
    task_id="Evaluate_biological_reproducibility", 
    bash_command='(cd %s; make {makeadd} Evaluate_biological_reproducibility) ' % dir, 
    dag=dag)

Evaluate_biological_reproducibility.set_upstream(Filter_lowly_expressed_Unigenes)
Evaluate_biological_reproducibility.set_downstream(Doc_assemble_quality)

Tabularize_Unigene_expression_annotation = BashOperator(
    task_id="Tabularize_Unigene_expression_annotation", 
    bash_command='(cd %s; make {makeadd} Tabularize_Unigene_expression_annotation) ' % dir, 
    dag=dag)

Tabularize_Unigene_expression_annotation.set_upstream(
        [Evaluate_biological_reproducibility, Summarize_annotation])
Tabularize_Unigene_expression_annotation.set_downstream(Doc_assemble_quality)


Screen_DE_Unigenes = BashOperator(
    task_id="Screen_DE_Unigenes", 
    bash_command='(cd %s; make {makeadd} Screen_DE_Unigenes) ' % dir, 
    dag=dag)

Screen_DE_Unigenes.set_upstream(Evaluate_biological_reproducibility)

Annotate_DE_Unigenes = BashOperator(
    task_id="Annotate_DE_Unigenes", 
    bash_command='(cd %s; make {makeadd} Annotate_DE_Unigenes) ' % dir, 
    dag=dag)

Annotate_DE_Unigenes.set_upstream([Screen_DE_Unigenes, Summarize_annotation])

Doc_DE_Unigenes = BashOperator(
    task_id="Doc_DE_Unigenes",
    bash_command='(cd %s; make {makeadd} d_c_de_gene d_d_de_gene_cluster) ' % dir,
    dag=dag)

Doc_DE_Unigenes.set_upstream([Annotate_DE_Unigenes])
Doc_DE_Unigenes.set_downstream(Doc_all)

Screen_DE_Isoforms = BashOperator(
    task_id="Screen_DE_Isoforms", 
    bash_command='(cd %s; make {makeadd} Screen_DE_Isoforms) ' % dir, 
    dag=dag)

Screen_DE_Isoforms.set_upstream(Evaluate_biological_reproducibility)

Annotate_DE_Isoforms = BashOperator(
    task_id="Annotate_DE_Isoforms", 
    bash_command='(cd %s; make {makeadd} Annotate_DE_Isoforms) ' % dir, 
    dag=dag)

Annotate_DE_Isoforms.set_upstream([Screen_DE_Isoforms, Summarize_annotation])


Doc_DE_Isoforms = BashOperator(
    task_id="Doc_DE_Isoforms",
    bash_command='(cd %s; make {makeadd} c_c_de_isoform c_d_de_isoform_cluster) ' % dir,
    dag=dag)

Doc_DE_Isoforms.set_upstream(Annotate_DE_Isoforms)
Doc_DE_Isoforms.set_downstream(Doc_all)

Success_mail = EmailOperator(
    task_id="Success_mail",
    to=cc_mail, 
    cc=mail, 
    subject="%s Finished" % prefix, 
    html_content='Dear Sir/Mandam<br> Please check your results in <a href="http://210.74.4.67:11521/result/{prefix}_EHBio">http://210.74.4.67:11521/result/{prefix}_EHBio</a> (username: ehbio; password: {prefix}_EHBio) <br> EHBIO Gene Technology'.format(prefix=prefix), 
    dag=dag)

Success_mail.set_upstream(Doc_all)
'''.format(mail=ehbioMailLs, customerMail=customerMailLs, prefix=prefix, 
        work_dir=work_dir, sampleL=sampleLs, makeadd=makeadd, owner=owner, 
        count = count)

    output_fh.close()

#-----airflow_etl_rnaseq_denovo----------------
def airflow_etl_rnaseq_ref(prefix, customerMailLs, ehbioMailLs, work_dir, sampleLs, dryRun, owner):
    if dryRun:
        makeadd = "-n"
    else:
        makeadd = ""
    count = 1 
    output = work_dir + '/' + prefix + "_"+str(count)+'_' + owner +'_ref_rnaseq_airflow_dag.py'
    while os.path.exists(output):
        print >>sys.stderr, output+" exists! Please check first."
        count += 1
        output = work_dir + '/' + prefix + "_"+str(count)+'_' + owner +'_ref_rnaseq_airflow_dag.py'
        #sys.exit(1)

    output_fh = open(output,  'w')
    print >>output_fh, '''   
from airflow.models import DAG
from airflow.operators.bash_operator import BashOperator
from airflow.operators.email_operator import EmailOperator

from datetime import datetime, timedelta


one_min_ago = datetime.combine(datetime.today() - timedelta(minutes=1),
                                  datetime.min.time())
#now = datetime.now()

mail = [{mail}]
cc_mail = [{customerMail}]
prefix = '{prefix}'
DAG_id = "Ref_RNAseq_{count}_" + prefix + "_{owner}"

dir = "{work_dir}"
nameL = [{sampleL}]

default_args = {{
    'owner': '{owner}',         
    'depends_on_past': False, 
    'start_date': one_min_ago, 
    'email': mail,
    'email_on_failure': True, 
    'email_on_retry': True, 
    'retries': 1000, 
    'retry_delay': timedelta(minutes=300), 
}}

dag = DAG(DAG_id, default_args=default_args,
    schedule_interval="@once")

Doc_all = BashOperator(
    task_id="Doc_all",
    bash_command='(cd %s; make {makeadd} DOC2) ' % dir,
    dag=dag) 

Doc_sequencing_quality = BashOperator(
    task_id="Doc_sequencing_quality", 
    bash_command='(cd %s; make {makeadd} a_main_fastqc_doc) ' % dir, 
    retry_delay=timedelta(hours=30),
    retries=100, 
    dag=dag)

Doc_sequencing_quality.set_downstream([Doc_all])

star_load_memory = BashOperator(
    task_id="star_load_memory", 
    bash_command="(cd %s; make {makeadd} star_load_memory) " % dir,
    dag=dag
    )


star_remove_memory = BashOperator(
    task_id="star_remove_memory", 
    bash_command="(cd %s; make {makeadd} star_remove_memory) " % dir,
    dag=dag
    )

map_summary = BashOperator(
    task_id="map_summary", 
    bash_command="(cd %s; make {makeadd} map_summary) " % dir,
    dag=dag
    )

geneBody_coverage_summary = BashOperator(
    task_id="geneBody_coverage_summary", 
    bash_command="(cd %s; make {makeadd} geneBody_coverage_summary) " % dir,
    dag=dag
    )

read_distrib_summary = BashOperator(
    task_id="read_distrib_summary", 
    bash_command="(cd %s; make {makeadd} read_distrib_summary) " % dir,
    dag=dag
    )

RPKM_saturation_summary = BashOperator(
    task_id="RPKM_saturation_summary", 
    bash_command="(cd %s; make {makeadd} RPKM_saturation_summary) " % dir,
    dag=dag
    )

Count_matrix = BashOperator(
    task_id="Count_matrix", 
    bash_command="(cd %s; make {makeadd} Count_matrix) " % dir,
    dag=dag
    )

TPM_matrix = BashOperator(
    task_id="TPM_matrix", 
    bash_command="(cd %s; make {makeadd} TPM_matrix) " % dir,
    dag=dag
    )

stringtie_merge = BashOperator(
    task_id="Stringtie_merge", 
    bash_command="(cd %s; make {makeadd} stringtie_merge) " % dir,
    dag=dag
    )

map_distrib_summary_doc = BashOperator(
    task_id="map_distrib_summary_doc", 
    bash_command="(cd %s; make {makeadd} map_distrib_summary_doc) " % dir,
    dag=dag
    )
map_distrib_summary_doc.set_upstream([map_summary, read_distrib_summary, geneBody_coverage_summary, RPKM_saturation_summary])

for i in nameL:
    SRA_download = BashOperator(
        task_id="Download_SRA_data_"+i, 
        bash_command=\
        '(cd %s; make {makeadd} %s.SRA_download sample=%s) ' % (dir, i, i), 
        retry_delay=timedelta(minutes=10),
        retries=100, 
        priority_weight=1000,
        dag=dag)

    SRA_rename = BashOperator(
        task_id="Rename_SRA_data_"+i, 
        bash_command=\
        '(cd %s; make {makeadd} %s.SRA_rename sample=%s) ' % (dir, i, i), 
        priority_weight=990,
        dag=dag)

    SRA_rename.set_upstream(SRA_download)

    Statistics_sequenced_reads_bases = BashOperator(
        task_id="Statistics_sequenced_reads_bases_"+i, 
        bash_command=\
        '(cd %s; make {makeadd} %s.statistics_sequenced_reads_bases sample=%s) ' \
            % (dir, i, i), 
        priority_weight=980,
        dag=dag)

    Statistics_sequenced_reads_bases.set_upstream(SRA_rename)
    Statistics_sequenced_reads_bases.set_downstream(Doc_sequencing_quality)

    fastqc = BashOperator(
        task_id="Quality_estimation_for_samp_"+i, 
        bash_command='(cd %s; make {makeadd} %s.fastqc sample=%s) ' \
            % (dir, i, i), 
        priority_weight=960,
        dag=dag)

    fastqc.set_upstream(SRA_rename)

    Trim = BashOperator(
        task_id="Trim_adaptors_lowquality_bases_for_samp_"+i, 
        bash_command='(cd %s; make {makeadd} %s.trim sample=%s) ' \
            % (dir, i, i), 
        priority_weight=940,
        dag=dag)

    Trim.set_upstream(fastqc)
    
    Split = BashOperator(
        task_id="Reads_equal_length_for_samp_"+i, 
        bash_command='(cd %s; make {makeadd} %s.split sample=%s) ' \
            % (dir, i, i), 
        priority_weight=920,
        dag=dag)

    Split.set_upstream(Trim) 
    Split.set_downstream([Doc_sequencing_quality])

    map = BashOperator(
        task_id="Reads_mapping_for_samp_"+i, 
        bash_command='(cd %s; make {makeadd} %s.map sample=%s) ' \
            % (dir, i, i), 
        priority_weight=900,
        dag=dag)
    map.set_upstream([Split, star_load_memory])
    map.set_downstream([star_remove_memory, map_summary, Count_matrix])

    stringtie_first = BashOperator(
        task_id="Transcript_assembling_for_samp_"+i, 
        bash_command='(cd %s; make {makeadd} %s.stringtie_first sample=%s) ' \
            % (dir, i, i), 
        priority_weight=900,
        dag=dag)
    stringtie_first.set_upstream([map])
    stringtie_first.set_downstream([stringtie_merge])

    geneBody_coverage = BashOperator(
        task_id="GeneBodyCoverage_estimating_for_samp_"+i, 
        bash_command='(cd %s; make {makeadd} %s.geneBody_coverage sample=%s) ' \
            % (dir, i, i), 
        priority_weight=880,
        dag=dag)
    geneBody_coverage.set_upstream([map])
    geneBody_coverage.set_downstream([geneBody_coverage_summary])

    read_distrib = BashOperator(
        task_id="Read_distrib_estimating_for_samp_"+i, 
        bash_command='(cd %s; make {makeadd} %s.read_distrib sample=%s) ' \
            % (dir, i, i), 
        priority_weight=880,
        dag=dag)
    read_distrib.set_upstream([map])
    read_distrib.set_downstream([read_distrib_summary])

    RPKM_saturation = BashOperator(
        task_id="RPKM_saturation_estimating_for_samp_"+i, 
        bash_command='(cd %s; make {makeadd} %s.RPKM_saturation sample=%s) ' \
            % (dir, i, i), 
        priority_weight=880,
        dag=dag)
    RPKM_saturation.set_upstream([map])
    RPKM_saturation.set_downstream([RPKM_saturation_summary])

    tpm = BashOperator(
        task_id="TPM_estimating_for_samp_"+i, 
        bash_command='(cd %s; make {makeadd} %s.tpm sample=%s) ' \
            % (dir, i, i), 
        priority_weight=880,
        dag=dag)
    tpm.set_upstream([map])
    tpm.set_downstream([TPM_matrix])


DE_genes = BashOperator (
    task_id="DE_genes", 
    bash_command='(cd %s; make {makeadd} DE_genes) ' % (dir), 
    retry_delay=timedelta(hours=240),
        dag=dag)
DE_genes.set_upstream([Count_matrix])


DE_genes_anno = BashOperator (
    task_id="DE_genes_anno", 
    bash_command='(cd %s; make {makeadd} DE_genes.anno) ' % (dir), 
    retry_delay=timedelta(hours=240),
        dag=dag)
DE_genes_anno.set_upstream([DE_genes])

DE_genes_go = BashOperator (
    task_id="DE_genes_go", 
    bash_command='(cd %s; make {makeadd} DE_genes.go) ' % (dir), 
    retry_delay=timedelta(hours=240),
        dag=dag)
DE_genes_go.set_upstream([DE_genes])

DE_genes_doc = BashOperator (
    task_id="DE_genes_doc", 
    bash_command='(cd %s; make {makeadd} DE_genes.doc) ' % (dir), 
    retry_delay=timedelta(hours=240),
        dag=dag)
DE_genes_doc.set_upstream([DE_genes, DE_genes_go, DE_genes_anno])
DE_genes_doc.set_downstream(Doc_all)


#Full_transcripts_annotation = BashOperator (
#    task_id="Full_transcripts_annotation", 
#    bash_command='(cd %s; make {makeadd} cuffcompare.anno) ' % (dir), 
#    retry_delay=timedelta(minutes=1),
#    dag=dag)
#
#Full_transcripts_annotation.set_upstream(Full_transcripts)
#
##DE_gene_cuffdiff = BashOperator (
##    task_id="DE_gene_cuffdiff", 
##    bash_command='(cd %s; make {makeadd} cuffdiff_cuffcompare) ' % (dir), 
##    dag=dag)
#
#DE_gene_DESeq2 = BashOperator (
#    task_id="DE_gene_DESeq2", 
#    bash_command='(cd %s; make {makeadd} deseq_cuffcompare.de) ' % (dir), 
#    retries=1000, 
#    retry_delay=timedelta(minutes=1),
#    dag=dag)
#
#for i in nameL:
#    Reads_cnt_htseq = BashOperator(
#        task_id="Reads_count_for_samp_"+i, 
#        bash_command='(cd %s; make {makeadd} %s.deseq_cuffcompare sample=%s) ' \
#            % (dir, i, i), 
#        retries=1000, 
#        retry_delay=timedelta(minutes=1),
#        priority_weight=840,
#        dag=dag)
#    Reads_cnt_htseq.set_upstream(Full_transcripts)
#    Reads_cnt_htseq.set_downstream(DE_gene_DESeq2)
#
#DE_gene_DESeq2_anno = BashOperator (
#    task_id="DE_gene_DESeq2_anno", 
#    bash_command='(cd %s; make {makeadd} deseq_cuffcompare.anno) ' % (dir), 
#    dag=dag)
#
#DE_gene_DESeq2_anno.set_upstream([DE_gene_DESeq2, Full_transcripts_annotation])
#
#DE_gene_enrichment = BashOperator (
#    task_id="DE_gene_enrichment", 
#    bash_command='(cd %s; make {makeadd} deseq_cuffcompare.go) ' % (dir), 
#    dag=dag)
#
#DE_gene_enrichment.set_upstream(DE_gene_DESeq2_anno)
#
#Doc_DE_gene = BashOperator (
#    task_id="Doc_DE_gene", 
#    bash_command='(cd %s; make {makeadd} c_de_gene_cuffcompare_doc) ' % (dir), 
#    dag=dag)
#
#Doc_DE_gene.set_upstream(DE_gene_enrichment)
#Doc_DE_gene.set_downstream(Doc_all)

#MATS = BashOperator (
#    task_id="Alternative_splicing_identification", 
#    bash_command='(cd %s; make {makeadd} MATS) ' % (dir), 
#    dag=dag)
#
#MATS.set_upstream(Full_transcripts)
#
#MATS_analyze = BashOperator (
#    task_id="Alternative_splicing_annotation", 
#    bash_command='(cd %s; make {makeadd} MATS_analyze) ' % (dir), 
#    dag=dag)
#
#MATS_analyze.set_upstream(MATS)
#MATS_analyze.set_upstream(Full_transcripts_annotation)
#
#Doc_AlternativeSplicing = BashOperator(
#    task_id="Doc_AlternativeSplicing",
#    bash_command='(cd %s; make {makeadd} h_mats_as) ' % dir,
#    dag=dag)
#
#Doc_AlternativeSplicing.set_upstream(MATS_analyze)
#Doc_AlternativeSplicing.set_downstream(Doc_all)

Success_mail = EmailOperator(
    task_id="Success_mail", 
    to=cc_mail,
    cc=mail, 
    subject="%s Finished" % prefix,  
    html_content='Dear Sir/Mandam<br> Please check your results in <a href="http://210.74.4.67:11521/result/{prefix}_EHBio">http://210.74.4.67:11521/result/{prefix}_EHBio</a> (username: ehbio; password: {prefix}_EHBio) <br> EHBIO Gene Technology'.format(prefix=prefix), 
    dag=dag)

Success_mail.set_upstream(Doc_all)

'''.format(mail=ehbioMailLs, customerMail=customerMailLs, prefix=prefix, 
        work_dir=work_dir, sampleL=sampleLs, makeadd=makeadd, owner=owner, 
        count=count)


#-------airflow_etl_rnaseq_ref---------------------    

def airflow_etl(hts_type, prefix, customerMailLs, ehbioMailLs, work_dir, sampleLs, dryRun, owner, assembl_only):
    if hts_type == "RNA_seq:denovo":
        airflow_etl_rnaseq_denovo(prefix, customerMailLs, ehbioMailLs, work_dir, sampleLs, dryRun, owner, assembl_only)
    elif hts_type == "RNA_seq:reference_based":
        airflow_etl_rnaseq_ref(prefix, customerMailLs, ehbioMailLs, work_dir, sampleLs, dryRun, owner)

#----airflow_etl-------------------------------

def checkAndUpdateMakefile(ppPath, hts_type, work_dir):
    fileD = {'RNA_seq:denovo':'Makefile.rna.noref.NEW2', 
             'RNA_seq:reference_based': 'Makefile.rna_std.NEW'}
    include = fileD.get(hts_type, '')
    if not include:
        print >>sys.stderr, "Pipeline for "+hts_type+' not included/'
        return
    addIn = 1
    if os.path.exists(work_dir+"/Makefile"):
        for line in open(work_dir+"/Makefile"):
            if line.find(include) != -1:
                addIn = 0
    #-------------------------------
    if addIn:
        fh_make = open(work_dir+"/Makefile", "a")
        print >>fh_make, "include "+ppPath.rstrip('/')+'/'+include
        fh_make.close()
#================================================

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    global makefileam

    file = options.filein
    hts_type = options.hts_type
    process = options.process
    AS = options.AS
    verbose = options.verbose
    assembl_only = options.assembl_only
    global debug
    debug = options.debug
    hts_typeL = options.hts_type.split(':')
    dryRun = options.dryRun
    ppPath = options.ppPath
    owner = options.owner
    #work_dir = options.work_dir
    work_dir = ''
    if not work_dir:
        work_dir = os.path.dirname(os.path.abspath(file))

    checkAndUpdateMakefile(ppPath, hts_type, work_dir)

    if options.check and os.path.exists(work_dir+'/makefile.am.template'):
        print >>sys.stderr, """
*********************************Warning******************************
*  File <makefile.am.template> exists. 
*  Please remove it and rerun the program to avoid losing information.
**********************************************************************

"""
        sys.exit(1)

    makefileam = open(work_dir+'/makefile.am.template', 'w')
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------

    context = ''.join(fh.readlines())
    contextD = json_loads(context)
    #print contextD.keys()

    #------------------------------------------------
    key = "Basic_parameter"
    traverseDict(key, contextD[key], 1, {}, hts_typeL)
    prefix = contextD[key]['prefix']
    checkLegalWord(prefix)

    customerMailL  = [i.strip() for i in contextD[key]['customerMail'].split(',')]
    #print >>sys.stderr, customerMailL
    ehbioMailL  = [i.strip() for i in contextD[key]['ehbioMail'].split(',')]
    if not customerMailL[0] or not ehbioMailL[0]:
        print >>sys.stderr, "Please define our mail address as well as customer mail address in config.json file."
        sys.exit(1)
    customerMailLs = grenerateQuotedLists(customerMailL)
    ehbioMailLs = grenerateQuotedLists(ehbioMailL)

    sampleD = contextD["Sample_info"]

    SRA_list = sampleD.get("SRA_list", '')
    if SRA_list:
        print >>makefileam, "\n#SRA information"
        for itemD in SRA_list:
            for key, valueL in itemD.items():
                srr = ' '.join(valueL)
                if srr.find('EHBIO') != -1:
                    print >>sys.stderr, "Wrong SRA information. \
Please delete `SRA_list` part and rerun the program. !!!"
                    sys.exit(1)
                print >>makefileam, "%s_sra=%s" % (key, srr)
        if hts_type.startswith("RNA_seq:denovo"):
            print >>makefileam, "\nadd_left_right_suffix=yes"
        print

    sampleL, sample_list = parseSampleInfo(sampleD, work_dir)
    sampleLs = grenerateQuotedLists(sampleL)

    sampleL.append(prefix)
    sampleCustomizeD = dict([(sample, {}) for sample in sampleL])

    #-------------------------------------------------------------------
    #for key, valueD in contextD.items():
    #    if key not in ignoreKeyD:
    #        level = 1
    #        traverseDict(key, valueD, level, sampleCustomizeD,
    #            hts_typeL)
    #-------------------------------------------------------------------
    if AS:
        print >>makefileam, "\n# If MATS==yes, reads in all samples will be trimmed to have similar length.\n"
        print >>makefileam, "\nMATS=yes\n"

    for key in ["HTS_common_parameter", "Analysis_type"]:
        traverseDict(key, contextD[key], 1, sampleCustomizeD, hts_typeL)
    
    #-------------Sample_specific_parameter---------
    Sample_specific_parameterD = contextD["Sample_specific_parameter"]
    for samp, valueD in Sample_specific_parameterD.items():
        if samp != "comment":
            sampD = sampleCustomizeD.get(samp, '')
            if sampD:
                for key, value in valueD.items():
                    if key == 'libtype':
                        parseLibType(sampD, value)
                    sampD[key] = value       
            elif samp.find('EHBIO') == -1:
                print >>sys.stderr, "Unexisting %s; \
                    Please check the spelling" % samp
                sys.exit(1)
    #-------------END reading file----------
    specialNameL = ['execute_htseq']
    specialNameD = {}
    print >>makefileam, "\n# Sample_specific_parameter\n"
    for samp in sampleL:
        sampD = sampleCustomizeD[samp]
        for key, value in sampD.items():
            print >>makefileam, "%s_%s=%s" % (samp, key, value)
            if samp != prefix and key in specialNameL:
                key = key.split("_",1)[1]
                if value != "no":
                    if key not in specialNameD:
                        specialNameD[key] = []
                    specialNameD[key].append(samp)
    #-------------------------------------------------
    print
    for key, itemL in specialNameD.items():
        print >>makefileam, "%s_nameL=%s" % (key, ' '.join(itemL))
    #----close file handle for files-----
    htseqL = specialNameD.get('htseq', [])
    sampleFileGeneration(sample_list, work_dir, htseqL)
    makefileam.close()
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    # Generate airflow etl

    airflow_etl(hts_type, prefix, customerMailLs, ehbioMailLs, work_dir, sampleLs, dryRun, owner, assembl_only)
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


