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
    This is designed to generate GO or other annotation file for goEnrichment.py.


Input:
  * Normally, at least 3 columns are needed, <Genes, PathwayID, PathwayDescrip>
  * Specially the program will work on only two columns <Genes, PathwayDescript>


Genes	KO Name	KO Description	Gene Ontology Biological Pathway	BP Description	KOG ID	KOG Description
c131183_g1	--	--	--	--	--	--
c69358_g1	4CL	4-coumarate--CoA ligase 	GO:0008152	metabolic process	KOG1176	Acyl-CoA synthetase
c142655_g1	"VCP, CDC48"	transitional endoplasmic reticulum ATPase	GO:0015979//GO:0006810//GO:0019079//GO:0015995//GO:0006281//GO:0015994//GO:0006310//GO:0070526//GO:0006576//GO:0006355	"photosynthesis//transport//viral genome replication//chlorophyll biosynthetic process//DNA repair//chlorophyll metabolic process//DNA recombination//threonylcarbamoyladenosine biosynthetic process//cellular biogenic amine metabolic process//regulation of transcription, DNA-templated"	KOG0730	AAA+-type ATPase


'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
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
        metavar="FILEIN", help="A matrix file as listed above")
    parser.add_option("-g", "--gene-name-column", dest="gene_col",
        default=1, type="int", 
        help="1-based number to spcify the column containing gene lists.")
    parser.add_option("-a", "--attribute-column", dest="attr_col",
        help="A list of attributes given as below <name1,id_col_num,descrip_col_num;name2,id_col_num,descript_col_num>. When using above example, we may give <BP,4,5;KEGG,6,7> to specify the column containing BP_id, BP_descrip, KEGG_id, KEGG_descrip. Specially <BP,4,4> will use BP_descrip as BP_id if BP_id not available.")
    parser.add_option("-s", "--separtor", dest="attrib_sep",
        default="//", help="As showed for gene <c142655_g1> in above example, it contains multiple GO annotations spearted by '//'. <'//'> should be given here as a separtor which is the default.")
    parser.add_option("-o", "--output-prefix", dest="output_p",
        help="Prefix for output file.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    gene_col = options.gene_col-1
    attr_colL = [i.split(',') for i in options.attr_col.split(';')]
    for attr_col in attr_colL:
        attr_col[1] = int(attr_col[1])-1
        attr_col[2] = int(attr_col[2])-1       
        attr_col.append({})
        attr_col.append(set())
        '''
        attr_col = [BP, id, descrip, 
                        {'id\tdescript':[gene1, gene2]}, set(gene1,gene2)}]
        '''

    attrib_sep = options.attrib_sep
    output_p   = options.output_p
    header = 1
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    for line in fh:
        if header:
            header -= 1
            continue
        lineL = line.strip().split('\t')
        gene = lineL[gene_col]
        for attr_col in attr_colL:
            id = lineL[attr_col[1]]
            desp = lineL[attr_col[2]].strip('"')
            if id == '--' and desp == '--':
                continue

            attr_col[4].add(gene)

            idL = id.split(attrib_sep)
            despL = desp.split(attrib_sep)
            len_idL = len(idL)

            for i in range(len_idL):
                id = idL[i]
                desp = despL[i]
                key = '\t'.join([id, desp])
                if key not in attr_col[3]:
                    attr_col[3][key] = [gene]
                else:
                    attr_col[3][key].append(gene)
    #-------------END reading file----------
    for attr_col in attr_colL:
        name = attr_col[0]
        output = output_p + '.' + name + '.xls'
        fh = open(output, 'w')
        total = str(len(attr_col[4]))
        print >>fh, "%s_id\t%s_description\tGene_list\tNo_gene_under_term\tTotal_annotated_gene" % (name, name)
        for key, valueL in attr_col[3].items():
            print >>fh, "\t".join([key, ','.join(valueL),str(len(valueL)), total])
        fh.close()
    #----close file handle for files-----
    if file != '-':
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
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


