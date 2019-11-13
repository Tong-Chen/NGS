#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import unicode_literals
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

tsv: 

1. Without header line

filename: data/00086b37-ad3a-4e4b-b44d-ea0cc657f48b/bdc49eae-31d4-425b-b7d1-49f1cf14df44.FPKM-UQ.txt.gz

(FPKM-UQ in filename will be used to label data)

ENSG00000000003.14  2129.48 
ENSG00000000004.14  1838.37 
ENSG00000000013.14  849.00 
ENSG00000000023.14  0.96

ENS2SYN

ENSG00000223972.5       DDX11L1
ENSG00000227232.5       WASH7P
ENSG00000278267.1       MIR6859-1
ENSG00000243485.3       RP11-34P13.3
ENSG00000274890.1       MIR1302-2
ENSG00000237613.2       FAM138A
ENSG00000268020.3       OR4G4P
ENSG00000240361.1       OR4G11P
ENSG00000186092.4       OR4F5


Metadata (48 columns file containing following infos)

file_name will be used to filter META data
No index

cases_0_submitter_id	cases_0_case_id	data_type	cases_0_samples_0_sample_type	cases_0_samples_0_tissue_type	file_name	cases_0_samples_0_submitter_id	cases_0_project_project_id	cases_0_samples_0_portions_0_analytes_0_aliquots_0_aliquot_id	cases_0_samples_0_sample_id	file_id	data_category	cases_0_samples_0_tumor_descriptor	cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id	analysis_workflow_type
TCGA-B3-4104	e15a108e-0b32-426f-926a-3fb1af0a82d3	Gene Expression Quantification	Primary Tumor		5079781b-14aa-4ff7-8ab8-b251e34a8eba.htseq.counts.gz	TCGA-B3-4104-01A	TCGA-KIRP	176143e0-adb3-4791-8075-2cbc16601cd5	bb0efb98-db70-499d-bc72-56af138c419c	fbd06307-dac8-460d-9fc5-2db29b43b35d	Transcriptome Profiling		TCGA-B3-4104-01A-02R-1351-07	HTSeq - Counts
TCGA-LN-A4A8	aa364b98-0097-43a4-976d-742faa0276df	Gene Expression Quantification	Primary Tumor		fbd32fa7-a9c7-4148-a3d4-69179c833f97.FPKM-UQ.txt.gz	TCGA-LN-A4A8-01A	TCGA-ESCA	1db4ac11-18c1-45ef-90de-fcb9d799d9ec	0c0582ca-a181-4efb-8486-c4312b250ba9	737a1837-5d14-45d7-af77-9b60fed2e556	Transcriptome Profiling		TCGA-LN-A4A8-01A-32R-A260-31	HTSeq - FPKM-UQ


1   caes_0_submitter_id
2   cases_0_case_id
3   data_type
4   cases_0_samples_0_sample_type
5   cases_0_samples_0_tissue_type
6   file_name
7   cases_0_samples_0_submitter_id
8   cases_0_project_project_id
9   cases_0_samples_0_portions_0_analytes_0_aliquots_0_aliquot_id
10  cases_0_samples_0_sample_id
11  file_id
12  data_category
13  cases_0_samples_0_tumor_descriptor
14  cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id
15  analysis_workflow_type

1   TCGA-B3-4104
2   e15a108e-0b32-426f-926a-3fb1af0a82d3
3   Gene Expression Quantification
4   Primary Tumor
5   
6   5079781b-14aa-4ff7-8ab8-b251e34a8eba.htseq.counts.gz
7   TCGA-B3-4104-01A
8   TCGA-KIRP
9   176143e0-adb3-4791-8075-2cbc16601cd5
10  bb0efb98-db70-499d-bc72-56af138c419c
11  fbd06307-dac8-460d-9fc5-2db29b43b35d
12  Transcriptome Profiling
13  
14  TCGA-B3-4104-01A-02R-1351-07
15  HTSeq - Counts

TCGA_barcode_sample_info

Description	cases_0_project_project_id	Cases_BCR	Cases_data	LastUpdate
Acute Myeloid Leukemia	LAML	200	200	05/31/16
Adrenocortical carcinoma	ACC	80	80	05/31/16
Bladder Urothelial Carcinoma	BLCA	412	412	05/27/16
Brain Lower Grade Glioma	LGG	516	516	05/02/16
Breast invasive carcinoma	BRCA	1100	1097	05/31/16
Cervical squamous cell carcinoma and endocervical adenocarcinoma	CESC	308	307	05/26/16
Cholangiocarcinoma	CHOL	36	36	05/31/16
Colon adenocarcinoma	COAD	461	461	05/27/16
Esophageal carcinoma	ESCA	185	185	05/31/16

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import h5py
import pandas as pd
from glob import glob

#from multiprocessing.dummy import Pool as ThreadPool

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
        metavar="FILEIN", help="A file folder containing TCGA tsv data. \
All files matching given patterns in this folder and subfolder will be used. \
Subfolder name should be <UUID>.")
    parser.add_option("-g", "--gene-list", dest="gene_list",
        help="A file contains list of gene symbols.")
    parser.add_option("-p", "--name-pat", dest="name_pat",
        default="*.txt.gz", help="Bash supported regular expression. Default <*.txt.gz>.")
    parser.add_option("-m", "--meta-data", dest="meta",
        help="TCGA Metadata file.")
    parser.add_option("-b", "--bar-code-samp", dest="barcode",
        default="/MPATHB/resource/TCGA/TCGA_barcode_sample_info.tsv", 
        help="TCGA sort code to full name description. Copied from TCGA website. Default </MPATHB/resource/TCGA/TCGA_barcode_sample_info.tsv>.")
    parser.add_option("-s", "--ensemble2genesymbol", dest="ens2syn",
        default="/MPATHB/resource/TCGA/anno/gencode.v22.ENS2SYN",  
        help="ENSEMBLE2GeneSymbol Map file. Default </MPATHB/resource/TCGA/anno/gencode.v22.ENS2SYN>.")
    parser.add_option("-o", "--output-prefix", dest="outp",
        help="Output prefix")
    parser.add_option("-t", "--save-type", dest="save_type",
        default='hdf5', help="TSV or hdf5 [default, fast, small, R-compatiable].")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    
    return (options, args)
#--------------------------------------------------------------------

def readTSV(tsv, exprD, gene_idL=[]):
    '''
    exprD = {'FPKM-QUAT': {'file1':matrix, 'file2':matrix}, 
             'FPKM':{'file1':matrix, 'file2':matrix}}
    tsv: /MPATHB/resource/TCGA/data/00086b37-ad3a-4e4b-b44d-ea0cc657f48b/bdc49eae-31d4-425b-b7d1-49f1cf14df44.FPKM-UQ.txt.gz
        00086b37-ad3a-4e4b-b44d-ea0cc657f48b: UUID
        bdc49eae-31d4-425b-b7d1-49f1cf14df44.FPKM-UQ.txt.gz: filename
        FPKM-UQ: type
    gene_idL: ENsemble gene ids. If given, only data about these genes will be extracted.
    '''
    uuid = os.path.split(tsv)[0].split('/')[-1]
    name = os.path.split(tsv)[-1]
    type = name.split('.')[1].replace('-', '_')
    #header=0: means treat the first line as header line
    ## data = pd.read_table(tsv, sep="\t", header=0)
    data = pd.read_table(tsv, sep="\t", index_col=0)
    # Set columns name
    data.columns = [uuid]

    # Set name of index and columns
    data.index.name = "gene_id"
    data.columns.name = 'uuid'

    # IN writing this func, <gene_start> is set to 'ENSG' in mind.
    # Filter rows in which <gene_id> column startws with <gene_start> 
    ## data = data[data['gene_id'].str.contains(r'^'+gene_start)]
    #data = data[data.index.str.contains(r'^'+gene_start)]
    if gene_idL:
        data = data[data.index.isin(gene_idL)]

    if type not in exprD:
        exprD[type] = {}
    exprD[type][name] = data
    #Select specific type
    #for type in typeL:
    #    assert name not in exprD[type], "Duplicate {}".format(name)
        ## expr_value = data[[gene_id, type]]
        ## use loc to avoid turning dataFrame to series
   #     expr_value = data.loc[:, [type]]
   #     expr_value.rename(index=str, columns={type:name}, inplace=True)
        
#-------------END readTSV---------------


def readMetaData(metadata):
    data = pd.read_table(metadata, sep="\t", header=0)
    #--get rownames and replace blank by '_'
    #colNames = [(name, name.replace(' ', '_')) for name in data.columns.values]
    #data.rename(columns=dict(colNames),  inplace=True)
    #return data

#-------------------------------------
def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    dir = options.filein
    pat  = options.name_pat
    meta_f = options.meta
    barcode_f = options.barcode
    outp = options.outp
    hdf5 = outp + '.hdf5'
    if options.gene_list:
        geneSymbolL = [i.strip() for i in open(options.gene_list)]
    else:
        geneSymbolL = []
    ens2syn = options.ens2syn
    verbose = options.verbose
    global debug
    debug = options.debug
    save_type = options.save_type
    #-----------------------------------
    ens2synDF = pd.read_table(ens2syn, sep="\t", header=0, index_col=0)
    if geneSymbolL:
        gene_idL = list(ens2synDF[ens2synDF["gene_symbol"].isin(geneSymbolL)].index)
        if debug:
            print >>sys.stderr, gene_idL
    else:
        gene_idL = []
    if save_type == 'hdf5':
        key = "ens2syn"
        hdf5_keys = [key]
        store = pd.HDFStore(hdf5, 'w', complib=str("zlib"), complevel=9)
        store[key] = ens2synDF
    #---------------------------------------------------------------
    
    # File in subdir of given directory
    fileL = []
    for dirpath, dirnames, filenames in os.walk(dir):
        fileL.extend(glob(os.path.join(dirpath, pat)))
            
    exprD = {}
    '''
    exprD = {'FPKM-UQ': {'file1':matrix, 'file2':matrix}, 
             'FPKM':{'file1':matrix, 'file2':matrix}}
    '''
    for file in fileL:
        if debug:
            print >>sys.stderr, file
        readTSV(file, exprD, gene_idL)

    meta = pd.read_table(meta_f, sep="\t", header=0)
    barcode = pd.read_table(barcode_f, sep="\t", header=0)

    typeL = exprD.keys()

    for type in typeL:
        filenameL = exprD[type].keys()
        exprMl = exprD[type].values()
        #exprMl.insert(ens2synDF)
        #Merge matrix by 'on' or by 'index'
        ## exprM_final = reduce(lambda left,right: pd.merge(left, right, on=gene_id, how='outer'), exprMl)
        exprM_final = reduce(lambda left,right: 
                pd.merge(left, right, left_index=True, right_index=True, how='outer'),
                exprMl)
        ## Substitute NaN to 0
        exprM_final = exprM_final.fillna(0)
        ## Remove all ZEROs
        exprM_final = exprM_final.loc[(exprM_final>0).any(axis=1)]

        ## exprM.to_csv(outp+'.'+type+'.xls', float_format="%.2f", index=False)
        ## b: bytes

        if save_type == 'hdf5':
            key = type
            hdf5_keys = [key]
            store[key] = exprM_final
            #exprM_final.to_hdf(hdf5, key, format="table", 
            #    mode='a', complevel=9, complib="bzip2")
        else:
            exprM_final.to_csv(outp+'.'+type+'.xls', sep=b"\t", float_format="%.2f")
        ## Add Gene symbol
        #exprM_withSym = pd.merge(ens2synDF, exprM_final, left_index=True, right_index=True, how="inner")
        #if save_type == 'hdf5':
        #    exprM_final.to_hdf(hdf5, type+'_withSym', format="table", 
        #        mode='a', complevel=9, complib="bzip2")
        #else:
        #    exprM_withSym.to_csv(outp+'.'+type+'.withSym.xls', sep=b"\t", 
        #        float_format="%.2f")
    
        #colN = exprM.columns.values()
        #newColN = meta.ix[colN]['Biosample term name']
    meta_part = meta[meta['file_name'].isin(filenameL)]
    meta_part = meta_part.fillna('')
    if save_type == 'hdf5':
        key = 'meta'
        hdf5_keys = [key]
        store[key] = meta_part
        store['barcode'] = barcode
        store.close()
        #meta_part.to_hdf(hdf5, key, format="table", 
        #    mode='a', complevel=9, complib="bzip2")
    else:
        meta_part.to_csv(outp+'.meta.xls', sep=b"\t")
        store.to_csv(outp+'.barcode.xls', sep=b"\t")
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
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


