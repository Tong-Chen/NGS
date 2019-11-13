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
    1. With header line (Normally the TPM and FPKM info will be extracted)
    2. Rows not start with <ENS> will be skipped.

gene_id transcript_id(s) length effective_length expected_count TPM FPKM posterior_mean_count posterior_standard_deviation_of_count pme_TPM pme_FPKM TPM_ci_lower_bound TPM_ci_upper_bound FPKM_ci_lower_bound FPKM_ci_upper_bound
8414    8414    93.00   0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0       0       0       0
8789    8789    82.00   0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0       0       0       0
ENSG00000000003.14 ENST00000373020.8,ENST00000494424.1,ENST00000496771.5,ENST00000612152.4,ENST00000614008.4 2129.48 1838.37 849.00 0.96 8.78 849.00  0.00    0.98    8.88    0.905918        1.04422 8.26776 9.52812

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

1   File accession
2   File format
3   Output type
4   Experiment accession
5   Assay
6   Biosample term id
7   Biosample term name
8   Biosample type
9   Biosample life stage
10  Biosample sex
11  Biosample organism
12  Biosample treatments
13  Biosample subcellular fraction term name
14  Biosample phase
15  Biosample synchronization stage
16  Experiment target
17  Antibody accession
18  Library made from
19  Library depleted in
20  Library extraction method
21  Library lysis method
22  Library crosslinking method
23  Experiment date released
24  Project
25  RBNS protein concentration
26  Library fragmentation method
27  Library size range
28  Biosample Age
29  Biological replicate(s)
30  Technical replicate
31  Read length
32  Mapped read length
33  Run type
34  Paired end
35  Paired with
36  Derived from
37  Size
38  Lab
39  md5sum
40  File download URL
41  Assembly
42  Platform
43  Controlled by
44  File Status
45  Audit WARNING
46  Audit INTERNAL_ACTION
47  Audit NOT_COMPLIANT
48  Audit ERROR

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
        metavar="FILEIN", help="A file folder containing ENCODE tsv data. All files matching given patterns in this folder will be used.")
    parser.add_option("-p", "--name-pat", dest="name_pat",
        default="ENC*.tsv", help="Bash supported regular expression Default <ENC*.tsv>.")
    parser.add_option("-m", "--meta-data", dest="meta",
        help="ENCODE Metadata file.")
    parser.add_option("-s", "--ensemble2genesymbol", dest="ens2syn",
        help="ENSEMBLE2GeneSymbol Map file")
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

def readTSV(tsv, exprD, gene_id='gene_id', gene_start='ENS', typeL=['TPM', 'FPKM']):
    '''
    exprD = {'TPM': {'file1':matrix, 'file2':matrix}, 
             'FPKM':{'file1':matrix, 'file2':matrix}}
    '''
    name = os.path.split(tsv)[-1][:-4]
    
    #header=0: means treat the first line as header line
    ## data = pd.read_table(tsv, sep="\t", header=0)
    data = pd.read_table(tsv, sep="\t", header=0, index_col=0)

    # IN writing this func, <gene_start> is set to 'ENSG' in mind.
    # Filter rows in which <gene_id> column startws with <gene_start> 
    ## data = data[data['gene_id'].str.contains(r'^'+gene_start)]
    data = data[data.index.str.contains(r'^'+gene_start)]
    
    #Select specific type
    for type in typeL:
        assert name not in exprD[type], "Duplicate {}".format(name)
        ## expr_value = data[[gene_id, type]]
        ## use loc to avoid turning dataFrame to series
        expr_value = data.loc[:, [type]]
        expr_value.rename(index=str, columns={type:name}, inplace=True)
        exprD[type][name] = expr_value
        
#-------------END readTSV---------------


def readMetaData(metadata):
    data = pd.read_table(metadata, sep="\t", header=0, index_col=0)
    #--get rownames and replace blank by '_'
    colNames = [(name, name.replace(' ', '_')) for name in data.columns.values]
    data.rename(columns=dict(colNames),  inplace=True)
    return data

#-------------------------------------
def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    dir = options.filein
    pat  = options.name_pat
    meta_f = options.meta
    outp = options.outp
    hdf5 = outp + '.hdf5'
    ens2syn = options.ens2syn
    verbose = options.verbose
    global debug
    debug = options.debug
    save_type = options.save_type
    #-----------------------------------
    gene_id='gene_id'
    gene_start='ENS'
    typeL=['TPM', 'FPKM']
    #-----------------------------------
    ens2synDF = pd.read_table(ens2syn, sep="\t", header=0, index_col=0)
    if save_type == 'hdf5':
        key = "ens2syn"
        hdf5_keys = [key]
        store = pd.HDFStore(hdf5, 'w', complib=str("zlib"), complevel=9)
        store[key] = ens2synDF
    fileL = glob(os.path.join(dir, pat))
    exprD = {}
    for type in typeL:
        exprD[type] = {}
    '''
    exprD = {'TPM': {'file1':matrix, 'file2':matrix}, 
             'FPKM':{'file1':matrix, 'file2':matrix}}
    '''
    for file in fileL:
        if debug:
            print >>sys.stderr, file
        readTSV(file, exprD, gene_id, gene_start, typeL)

    meta = readMetaData(meta_f)

    for type in typeL:
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
    #meta_part = meta[meta.index.isin(exprM_final.T.index)]
    meta_part = meta[meta.index.isin(exprM_final.columns.values)]
    if save_type == 'hdf5':
        key = 'meta'
        hdf5_keys = [key]
        store[key] = meta_part
        store.close()
        #meta_part.to_hdf(hdf5, key, format="table", 
        #    mode='a', complevel=9, complib="bzip2")
    else:
        meta_part.to_csv(outp+'meta.xls', sep=b"\t")

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


