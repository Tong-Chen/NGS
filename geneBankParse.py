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
    This is designed to extract protein and CDS sequences from gene bank file named with protein_id.

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
from Bio import SeqIO

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
        metavar="FILEIN", help="GeneBank file")
    parser.add_option("-d", "--id-prefix", dest="id_prefix", default='', 
        help="Supply a string to output sequence IDs. Default no extra prefix.")
    parser.add_option("-G", "--genome", dest="genome",
        default=False, action="store_true", 
        help="Extract Genome sequence and gene position file.")
    parser.add_option("-g", "--gene", dest="gene",
        default=False, action="store_true", 
        help="Extract Gene sequence")
    parser.add_option("-c", "--cds", dest="cds",
        default=False, action="store_true", help="Extract CDS sequence")
    parser.add_option("-p", "--prot", dest="prot",
        default=False, action="store_true", help="Extract protein sequence")
    parser.add_option("-I", "--include_organism", dest="include_organism",
        default=False, action="store_true", help="Include organism information in output IDs")
    parser.add_option("-o", "--output-prefix", dest="output_prefix",
        help="Prefix of output files.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    genome = options.genome
    gene = options.gene
    cds = options.cds
    prot = options.prot
    id_prefix = options.id_prefix
    include_organism = options.include_organism
    output_prefix = options.output_prefix
    if genome:
        genome_fh =open(output_prefix+'.genome.fa', 'w') 
        gene_pos_fh = open(output_prefix+'.gene_pos.bed', 'w')
    if gene:
        gene_fh = open(output_prefix+'.gene.fa', 'w')
    if cds:
        cds_fh = open(output_prefix+'.cds.fa', 'w')
    if prot:
        prot_fh = open(output_prefix+'.prot.fa', 'w')
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    strandD = {1:'+'}
    genebank_parse = SeqIO.parse(file, 'genbank')
    for each_seq in genebank_parse:
        each_seq_annoD = each_seq.annotations
        organism = each_seq_annoD['organism']
        contig_name = each_seq.name
        if genome:
            genome_seq = each_seq.seq 
            print >>genome_fh, ">{}\n{}".format(contig_name, genome_seq)
        if include_organism:
            id_suffix = '@'+organism+'@'+contig_name
        else:
            id_suffix = ''
        
        for each_feature in each_seq.features:
            if each_feature.type == 'gene':
                if 'gene' in each_feature.qualifiers:
                    gene_id = each_feature.qualifiers['gene'][0]
                else:
                    gene_id = each_feature.qualifiers['locus_tag'][0]
                if gene:
                    print >>gene_fh, ">{}{}{}".format(id_prefix, gene_id, id_suffix)
                    print >>gene_fh, each_feature.location.extract(each_seq).seq
                if genome:
                    start = each_feature.location.start.position
                    end = each_feature.location.end.position
                    strand = each_feature.location.strand
                    strand = strandD.get(strand, '-')
                    print >>gene_pos_fh, "\t".join([contig_name, str(start), str(end), gene_id, organism, strand])
            #---------------------------------------------------------

            if (cds or prot) and each_feature.type == "CDS":
                prot_id = each_feature.qualifiers['protein_id'][0]
                if cds:
                    print >>cds_fh, ">{}{}".format(id_prefix, prot_id)
                    print >>cds_fh, each_feature.location.extract(each_seq).seq
                if prot:
                    print >>prot_fh, ">{}{}".format(id_prefix, prot_id)
                    print >>prot_fh, each_feature.qualifiers['translation'][0]
            #----------END each_feature-------------------------------
        #----------END each_seq-------------------------------
    #------------END all-------------------
    if cds:
        cds_fh.close()
    if prot:
        prot_fh.close()
    if gene:
        gene_fh.close()
    if genome:
        genome_fh.close()
        gene_pos_fh.close()
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


