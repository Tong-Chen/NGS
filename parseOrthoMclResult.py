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
    This is designed to parse orthmcl results.

    Input file format:
    
    cluster_name<colon><any blank>spe1<vertical_line>prot1<any blank>spe2<verticial_line>prot2<any blank>.....

    C10000: Aco|Aco000153.1 Aco|Aco004369.1 Aco|Aco010005.1
    C10001: Aco|Aco000153.1 Cla|Cla004369.1 Dec|Dec010005.1

    Tasks:

    1. Get a matrix showing the number of proteins in each cluster.
    2. Extract single gene clusters and their sequences in all given
    species. In the output nucleotide file, ending stop codon (TAA,
    TAG, TGA) will be removed for compatible with
    `translatorx_vLocal.pl` and `trimal`.
    3. Extract species specific clusters for given species.
    4. Extract gene-expansion clusters for given species.
    5. Extract multiple-species specific clusters.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
import re

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
        metavar="FILEIN", help="Output of `orthomclMclToGroups`.")
    parser.add_option("-t", "--target-species", dest="main_spe",
        help="Specify the `species` name used for extracting \
species specific clusters or specially expanded clusters.")
    parser.add_option("-E", "--exclude-all",
        dest="exclude_when_reading",
        help="Comma or blank separated strings representing \
species excluded when reading in the result. It will affect \
all tasks. Default including all species.")
    parser.add_option("-e", "--exclude-2", dest="exclude_single_conserve",
        help="Comma or blank separated strings representing \
species should not be considered when performing task <2>. \
Default including all species.")
    parser.add_option("-s", "--specific-multiple-5", 
        dest="specific_multiple",
        help="Comma or blank separated strings representing \
multiple species used for task <5>. \
Default muting task 5.")
    parser.add_option("-P", "--directory-prot", dest="dir_prot",
        help="Directory containing all protein sequences used \
for `orthoMcl.sh`. All sequences have a suffix `.fasta`.")
    parser.add_option("-N", "--directory-nucl", dest="dir_nucl",
        help="Directory containing all nucleotide sequences used \
for `orthoMcl.sh`. All sequences have a suffix `.fasta`.")
    parser.add_option("-o", "--output-prefix", dest="outp",
        help="Prefix for output files.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def getSpecific(speciesD, speL):
    """ 
    speciesD = {spe_1: {cluster_1:1, cluster_2:3, ...}, 
                spe_2: {cluster_2:2, cluster_3:4, ...}
               }
    """
    main_spe = speL[0]
    target_cluster_set = set(speciesD[main_spe].keys())
    for tmp_spe in speL[1:]:
        target_cluster_set.intersection_update(set(speciesD[tmp_spe].keys()))
    #--------------------------------------
    spe_allL = speciesD.keys()
    for tmp_spe in spe_allL:
        if tmp_spe not in speL:
            target_cluster_set.difference_update(set(speciesD[tmp_spe].keys()))
    #------------------------------------
    return target_cluster_set
#--------------getSpecific---------------------

def readFastaDir(dir):
    fileL = os.listdir(dir)
    seqD = {}
    for file in fileL:
        if not file.endswith('fasta'):
            continue
        key = ''
        for line in open(dir+'/'+file):
            if line[0] == '>':
                if key:
                    seqD[key] = ''.join(seqL)
                key = line[1:-1]
                assert key not in seqD, "Duplicate key in %s" % dir
                seqL = []
            else:
                seqL.append(line.strip())
        #--------------------------------------
        if key:
            seqD[key] = ''.join(seqL)
        #--------------------------------------
    return seqD
#-------------readDir-----------------

def getClusterSeq(clusterD, protseqD, nuclD, excludeSpeciesL, outp):
    '''
    clusterD = {cluster_1: {spe_1:[id1, id2, ...], 
                            spe_2:{id1, id2, ...}}, 
                cluster_2: {}, 
               }
    '''
    outp = outp + '.conserved_single_member_orthlog/'
    cmd = "mkdir -p " + outp
    os.system(cmd)

    for cluster, memD in clusterD.items():
        if protseqD:
            fhP = open(outp+cluster+'.pep.fa', 'w')
        if nuclD:
            fhN = open(outp+cluster+'.nucl.fa', 'w')
        #----------------------------------------------
        for spe, memL in memD.items():
            if spe in excludeSpeciesL:
                continue
            #--------------------------------
            if protseqD:
                for prot in memL:
                    print >>fhP, ">%s_%s\n%s" % \
                        (cluster, prot, protseqD[prot])
            #-------------------------------------
            if nuclD:
                for nucl in memL:
                    nucl_seq = nuclD[nucl]
                    nucl_seq_3 = nucl_seq[-3:]
                    if nucl_seq_3 in ["TAA", "TAG", "TGA"]:
                        nucl_seq = nucl_seq[:-3]
                    print >>fhN, ">%s_%s\n%s" % \
                        (cluster, nucl, nucl_seq)
                    #exclude ending stop codon, TAA, TAG, TGA
            #-------------------------------------
        #----------------------------------------
        if protseqD: fhP.close()
        if nuclD: fhN.close()
    #-----------END for-----------------
#---------getClusterSeq-------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    main_spe = options.main_spe
    outp = options.outp + '.orthomcl'
    dir_prot = options.dir_prot
    dir_nucl = options.dir_nucl
    #-----------------------------------------
    exclude_when_reading = options.exclude_when_reading
    exclude_when_readingL = []
    if exclude_when_reading:
        exclude_when_readingL = re.split('[, ]*', exclude_when_reading)
    #-----------------------------------------
    exclude_single_conserve = options.exclude_single_conserve
    exclude_single_conserveL = []
    if exclude_single_conserve:
        exclude_single_conserveL = re.split('[, ]*', exclude_single_conserve)
    #---------------------------------------------
    specific_multiple = options.specific_multiple
    specific_multipleL = []
    if specific_multiple:
        specific_multipleL = re.split('[, ]*', specific_multiple)
    #-----------------------------------------
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    protD = {}
    nuclD = {}
    if dir_prot:
        protD = readFastaDir(dir_prot)
    if dir_nucl:
        nuclD = readFastaDir(dir_nucl)
    if debug:
        print >>sys.stderr, protD
    #--------------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    matrixD = {}
    """ 
    matrixD = {cluster_1: {spe_1:5, spe_2:4, ...}, 
               cluster_2: {spe_1:6, spe_2:3, ...},
              }
    """
    speciesD = {}
    """ 
    speciesD = {spe_1: {cluster_1:1, cluster_2:3, ...}, 
                spe_2: {cluster_2:2, cluster_3:4, ...}
               }
    """
    clusterD = {}
    '''
    clusterD = {cluster_1: {spe_1:[id1, id2, ...], 
                            spe_2:{id1, id2, ...}}, 
                cluster_2: {}, 
               }
    '''

    for line in fh:
        try:
            cluster, protein = line.split(':')
        except ValueError:
            print >>sys.stderr, line
            sys.exit(1)
        proteinL = protein.strip().split()
        assert cluster not in matrixD, "Duplicate cluster %s" % cluster
        matrixD[cluster] = {}
        clusterD[cluster] = {}
        tmpD = {}
        clusterTmpD = {}
        for protein in proteinL:
            spe, prot = protein.split('|')
            if spe in exclude_when_readingL:
                continue
            tmpD[spe] = tmpD.get(spe, 0) + 1
            if spe not in clusterTmpD:
                clusterTmpD[spe] = []
            clusterTmpD[spe].append(protein)
            if spe not in speciesD:
                speciesD[spe] = {}
            speciesD[spe][cluster] = speciesD[spe].get(cluster, 0)+1
        matrixD[cluster] = tmpD
        clusterD[cluster] = clusterTmpD
    #-------------END reading file----------
    if debug:
        print >>sys.stderr, matrixD
        print >>sys.stderr, speciesD
    #----------------------------------------------
    speL = speciesD.keys()
    speL.sort()
    #------output matrix and get single protein family conserve---------------
    matrix = outp + '.matrix.xls'

    #single_conserve_prot = outp + '.conserved_single_member_orthlog_prot.fa'
    #single_conserve_prot_fh = open(single_conserve_prot, 'w')

    #single_conserve_nucl_fh = ''
    #if nuclD:
    #    single_conserve_nucl = outp + '.conserved_single_member_orthlog_nucl.fa'
    #    single_conserve_nucl_fh = open(single_conserve_nucl, 'w')
    single_conserve_speL = [tmp_spe for tmp_spe in speL \
        if tmp_spe not in exclude_single_conserveL]
    single_conserve_clusterD = {}
    
    #print >>single_conserve_fh, \
    #    "#species list:%s" % '\t'.join(single_conserve_speL)

    matrix_fh = open(matrix, 'w')
    print >>matrix_fh, "cluster\t%s" % '\t'.join(speL)
    for cluster,tmpD in matrixD.items():
        tmpL = [str(tmpD.get(tmp_spe, 0)) for tmp_spe in speL]
        print >>matrix_fh, "%s\t%s" % (cluster, '\t'.join(tmpL))
        
        #-----get conserved single member family-------
        saveSingle = 1
        for tmp_spe in single_conserve_speL:
            if tmpD.get(tmp_spe,  0) != 1:
                saveSingle = 0
                break
        #-----saveSingle=0 exclude--------
        if saveSingle:
            single_conserve_clusterD[cluster] = clusterD[cluster]
        #-----get conserved single member family-------
    matrix_fh.close()

    #-----get conserved single member family-------
    getClusterSeq(single_conserve_clusterD, protD, nuclD, exclude_single_conserveL, \
        outp)
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    
    #-------main specific-----------------------------
    main_spe_specific_clusterL = getSpecific(speciesD, [main_spe])
    if debug:
        print >>sys.stderr, main_spe_specific_clusterL
    #------------Output----------------------------
    main_spe_specific_file = outp + '.' + main_spe+'.specific_cluster.xls'
    main_spe_specific_fh = open(main_spe_specific_file, 'w')
    for cluster in main_spe_specific_clusterL:
        idL = [i for j in clusterD[cluster].values() for i in j]
        print >>main_spe_specific_fh, "%s:%s" % (cluster, '\t'.join(idL))
    main_spe_specific_fh.close()

    #-------multiple specific-----------------------------
    if specific_multipleL:
        specific_multipleL_cluster = \
            getSpecific(speciesD,specific_multipleL)
        if debug:
            print >>sys.stderr, specific_multipleL_cluster
        #------------Output----------------------------
        specific_multiple_file = outp + '.' + \
            '_'.join(specific_multipleL)+'.specific_cluster.xls'
        specific_multiple_fh = open(specific_multiple_file, 'w')
        for cluster in specific_multipleL_cluster:
            idL = [i for j in clusterD[cluster].values() for i in j]
            print >>specific_multiple_fh, "%s:%s" % (cluster, '\t'.join(idL))
        specific_multiple_fh.close()
    #-------multiple specific-----------------------------
    


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


