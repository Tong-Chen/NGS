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
    Simplify Trinotate annotation and summary the annotation results.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from bs4 import BeautifulSoup

#reload(sys)
#sys.setdefaultencoding('utf8')

#from multiprocessing.dummy import Pool as ThreadPool

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
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="The output of \
Trinotate/Trinotate.sqlite report.")
    parser.add_option("-g", "--uniprot2gene2reactome", dest="uniprot2gene",
        default="/MPATHB/resource/Reactome/Uniprot_symbol_reactome_withdef.txt", 
        help="Default /MPATHB/resource/Reactome/Uniprot_symbol_reactome_withdef.txt")
    parser.add_option("-k", "--ko2pathway", dest="ko2pathway",
        default="/MPATHB/resource/KEGG/KEGG.ko.pathway.txt", 
        help="Default /MPATHB/resource/KEGG/KEGG.ko.pathway.txt")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def getUniprotId(annofile):
    blastL = ['sprot_Top_BLASTX_hit',
        'sprot_Top_BLASTP_hit',
        'TrEMBL_BLASTX', 'TrEMBL_BLASTP']
    header = 1
    uniprotIDs = set([])
    for line in open(annofile):
        lineL = line.strip().split('\t')
        if header:
            indexL = [lineL.index(i) for i in blastL if i in lineL]
            header -= 1
            continue
        #------------------------------------
        idL = set([lineL[i].split('^', 1)[0] for i in indexL if lineL[i]!='.'])
        uniprotIDs.update(idL)
    return uniprotIDs
#------------------------------------------

def process(item, label, tmpL, header, uniprot2geneD, ko2pathwayD):
    blastL = ['sprot_Top_BLASTX_hit', 'TrEMBL_Top_BLASTX_hit',
        'sprot_Top_BLASTP_hit', 'TrEMBL_Top_BLASTP_hit',
        'TrEMBL_BLASTX', 'TrEMBL_BLASTP']
    #'Pfam', 'SignalP', 'TmHMM', 'eggnog',
    eggnog = ['eggnog']
    pfam = ['Pfam']
    go = ['gene_ontology_blast', 'gene_ontology_pfam']
    kegg = ['Kegg']
    depleted = ['RNAMMER', 'prot_coords', 'transcript', 'peptide']
    directReturn = ['prot_id','eggnog','SignalP','TmHMM']
    if label in blastL:
        if header:
            tmpL.append(item+'_id')
            tmpL.append(item+'_desp')
        else:
            if item == '.':
                tmpL.append('.')
                tmpL.append('.')
            else:
                itemL = item.split('^')
                match_id = itemL[0]
                #if match_id != '.' and tmpL[2] == "NA":
                #    tmpL[2] = 
                tmpL.append(match_id)
                tmpL.append(itemL[5].replace('Full=',\
                    '').split('{')[0].\
                    replace('RecName: ', "").replace("SubName: ", ""))
    elif label in kegg:
        if header:
            tmpL.append("KEGG ID")
            tmpL.append("KO")
            tmpL.append("KEGG Pathway")
            tmpL.append("Reactome")
        else:
            tmpL.append(item)
            itemL = item.split('`')
            kegg_pathway = '.'
            ko = '.'
            if len(itemL) >= 2:
                ko = itemL[-1]
                kegg_pathway = '`'.join(ko2pathwayD.get(ko, ['.']))
            #elif len(itemL) > 2:
            #    print >>sys.stderr, "Unexpected KEGG", item
            #    sys.exit(1)
            tmpL.append(ko)
            tmpL.append(kegg_pathway)
            tmpL.append('.')
    elif label in go:
        tmpL.append(item)
    elif label in directReturn:
        tmpL.append(item)
    elif label in pfam:
        if header:
            tmpL.append(item)
        else:
            if item == '.':
                tmpL.append('.')
            else:
                itemL = item.split('`')
                valueTl = []
                for item in itemL:
                    pfID, pfshort, pflong = item.split('^')[:3]
                    pfID = pfID.rsplit('.', 1)[0]
                    valueTl.append('^'.join([pfID, pfshort, pflong]))
                tmpL.append('`'.join(valueTl))
                #tmpL.append(itemL[0])
                #tmpL.append(itemL[1])
                #tmpL.append(itemL[2])
    elif label in depleted:
        pass
#---------------------------------
def save(annoD, labelD, gene, item, label):
    '''
    
    '''
    if item == '.': return

#--------------------------------------

def readko2pathway(ko2pathway):
    '''
    PathwayName     koName  PathwayDef      KOgene secKEGG firstKEGG
    '''

    ko2pathwayD = {}
    header = 1
    for line in open(ko2pathway):
        pathid, koid, pathdef, koname, secKEGG, firstKEGG = line.strip().split('\t')
        pathid = pathid.replace('path:', '')
        koid = koid.upper()
        if koid not in ko2pathwayD:
            ko2pathwayD[koid] = []
        ko2pathwayD[koid].append('^'.join([pathid, pathdef, secKEGG, firstKEGG]))
    return ko2pathwayD
#-------------------------------------------------

def readuniprot2gene(uniprot2gene, uniprotIDs):
    uniprot2geneD = {}
    header = 1
    for line in open(uniprot2gene):
        if header:
            header -= 1
            continue
        #---------------------
        name, UniProtKBID, Gene_Name, Reactome = line.strip().split('\t')
        if uniprotIDs and (UniProtKBID not in uniprotIDs):
            continue
        if UniProtKBID not in uniprot2geneD:
            uniprot2geneD[UniProtKBID] = {}
            uniprot2geneD[UniProtKBID]['Gene_Name'] = Gene_Name
            uniprot2geneD[UniProtKBID]['Reactome']  = Reactome
        else:
            print >>sys.stderr, "Duplicate", UniProtKBID
            sys.exit(1)
    #--------------------------------
    return uniprot2geneD
#------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    uniprot2gene = options.uniprot2gene
    uniprot2geneD = {}
    ko2pathway   = options.ko2pathway
    ko2pathwayD  = {}
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    annoD = {} #annoD = {"sprot":["g1", "g2"], "TrEMBL":["g1", "g3"]}
    labelD = {"sprot_Top_BLASTX_hit": "SwissProt",
            "TrEMBL_Top_BLASTX_hit":"TrEMBL", 
            "TrEMBL_BLASTX":"TrEMBL", 
            "TrEMBL_BLASTP":"TrEMBL", 
            "sprot_Top_BLASTP_hit":"SwissProt",
            "TrEMBL_Top_BLASTP_hit":"TrEMBL", "Pfam":"Pfam",
            "SignalP":"SignalP", "Kegg":"KEGG"}
    
    if ko2pathway:
        ko2pathwayD = readko2pathway(ko2pathway)
    #-------------------------------------------
    if uniprot2gene:
        uniprotIDs = set([])
        uniprotIDs = getUniprotId(file)
        uniprot2geneD = readuniprot2gene(uniprot2gene, uniprotIDs)

    header = 1
    headerL = []
    symbolD = {}
    for line in fh:
        if debug:
            print >>sys.stderr, line, 
        tmpL = []
        if header:
            headerL = line.strip().split('\t')
            tmpL.extend(headerL[:2])
            len_header = len(headerL)
            header -= 1
            for i in range(2, len_header):
                item  = headerL[i]
                label = headerL[i]
                process(item, label, tmpL, 1, uniprot2geneD, ko2pathwayD)
            tmpL.insert(2, "Gene Symbol")
            tmpL.insert(3, "Description")
            tmpL.insert(4, "Matched IDs")
            tmpL.append('GeneOntology')
            if debug:
                print >>sys.stderr, tmpL
            #headerL = tmpL[:]
            tmpL[0] = "Gene ID"
            tmpL[1] = "Transcript ID"
            tmpL[7] = "Protein ID"
            tmpL[14] = "PFAM"
            tmpL[17] = "EggNOG"
            #tmpL[18] = "KEGG ID"
            #tmpL[19] = "KEGG Description"
            tmpL[24] = "Gene Ontology"
            print '\t'.join(tmpL)
            if debug:
                print >>sys.stderr, headerL
            continue
        lineL = line.strip().split('\t')
        gene = lineL[0]
        annoD[gene] = set(['All'])
        #tr   = lineL[1]
        tmpL.extend(lineL[:2])
        # add gene name column
        tmpL.append('NA')
        tmpL.append('.')
        tmpL.append('.')
        otuputL = lineL[:2]
        for i in range(2, len_header):
            item  = lineL[i]
            label = headerL[i]
            process(item, label, tmpL, 0, uniprot2geneD, ko2pathwayD)
            if item != '.' and label in labelD:
                annoD[gene].add(labelD[label])
        tmpL.append('.')
        if debug:
            print >>sys.stderr, len(tmpL), tmpL
        #----get gene symbol and final description----
        # Columns for sprot_Top_BLASTX_hit_id, sprot_Top_BLASTX_hit_desp
        # sprot_Top_BLASTP_hit_id, sprot_Top_BLASTP_hit_desp
        # TrEMBL_BLASTX_id, TrEMBL_BLASTX_desp, TrEMBL_BLASTP_id, TrEMBL_BLASTP_desp
        id_despPair = [(8, 9), (5, 6), (12, 13), (10, 11)]
        for id_index, desp_index in id_despPair:
            id1 = tmpL[id_index]
            desp = tmpL[desp_index]
            if tmpL[2] == 'NA' and id1 != '.':
                symbol = uniprot2geneD.get(id1, {}).get("Gene_Name", 'NA')
                symbol = symbol.strip(';').upper()
                if symbol != 'NA' and desp != '.':
                    if symbol not in symbolD:
                        symbolD[symbol] = set([tmpL[0]])
                    else:
                        symbolD[symbol].add(tmpL[0])
                    symbol = symbol+'_'+str(len(symbolD.get(symbol)))
                    tmpL[2] = symbol
                    tmpL[3] = desp.replace(';', '')
                    tmpL[4] = id1
                    tmpL[21] = uniprot2geneD.get(id1, {}).get("Reactome", '.')
                    break
        #---------------------------------------------------
        if tmpL[2] == 'NA':
            for id_index, desp_index in id_despPair:
                id1 = tmpL[id_index]
                desp = tmpL[desp_index]
                if tmpL[2] == 'NA' and id1 != '.':
                    symbol = id1.split('_')[0]
                    if symbol != 'NA' and desp != '.':
                        if symbol not in symbolD:
                            symbolD[symbol] = set([tmpL[0]])
                        else:
                            symbolD[symbol].add(tmpL[0])
                        symbol = symbol+'_'+str(len(symbolD.get(symbol)))
                        tmpL[2] = symbol
                        tmpL[3] = desp.replace(';', '')
                        tmpL[4] = id1
                        tmpL[21] = uniprot2geneD.get(id1, {}).get("Reactome", '.')
                        break
        #-------------------------------------------------
        if tmpL[2] == 'NA':
            tmpL[2] = tmpL[0]
        if tmpL[21] == 'NA':
            tmpL[21] = '.'
        #---get gene symbol and final description----
        #---parse GO---------------------------------
        # Column for gene_ontology_blast, gene_ontology_pfam
        goL = []
        go = '.'
        go_blast = tmpL[-3]
        go_pfam  = tmpL[-2]
        if go_blast != '.':
            goL.extend([i for i in go_blast.split('`')])
        if go_pfam != '.':
            goL.extend([i for i in go_pfam.split('`')])
        if goL:
            goL = list(set(goL))
            goL = [i.split('^') for i in goL]
            goL = ['^'.join([i[1].capitalize().replace('_', ' '), i[0], i[2]]) for i in goL]
            goL.sort()
            go = '`'.join(goL)
        tmpL[-1] = go
        #---parse GO---------------------------------
        if debug:
            print >>sys.stderr, tmpL
        #--------------------------------------------
        print '\t'.join(tmpL)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    #--------Output statistis information------
    anno_fp = open(file+".sta.xls", 'w')
    for gene, labelL in annoD.items():
        for label in labelL:
            print >>anno_fp, "%s\t%s" % (gene, label)
    anno_fp.close()
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


