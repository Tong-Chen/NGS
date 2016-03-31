#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
import sys
from ncbi import taxonID2Name

def gi2taxon(gi_taxid_prot, taxon_id_name):
    aDict = {}
    for line in open(gi_taxid_prot):
        gi, taxid = line.strip().split()
        if taxid in taxon_id_name:
            taxname = taxon_id_name[taxid]
            aDict[gi] = taxname
        else:
            print >>sys.stderr, "No corresponding name %s" % taxid
            aDict[gi] = taxid
    #-------------------------------
    return aDict
#---------------------------------------


#def gitaxon(hitseq):
#    aDict = {}
#    for line in open(hitseq):
#        if line[0] == '>':
#            lineL = line.split('|',2)
#            gi = lineL[1]
#            taxonOri = lineL[2]
#            start = taxonOri.find('[')
#            #assert start != -1, taxonOri
#            end   = taxonOri.find(']')
#            #assert end != -1, taxonOri
#            if start != -1 and end != -1:
#                taxon = taxonOri[start+1:end]
#            else:
#                print >>sys.stderr, "%s Wrong taxon" % gi
#                taxon = 'NoTaxon'
#            aDict[gi] = taxon
#    #---------------------------------------
#    return aDict
##------------ENd gitaxon------------------------

def giToID(transferFile):
    aDict = {}
    for line in open(transferFile):
        gi, locus = line.strip().split()
        aDict[gi] = locus
    return aDict
#----------------------------------------

def main():
    print >>sys.stderr, "Print the result to screen"
    print >>sys.stderr, "Take the result from psiblastOneMsaOneTable,\
and add taxon information to discrimate the conservation of \
repetitions."
    if len(sys.argv) != 5:
        print >>sys.stderr, 'Using python %s msaresult hitseq \
giToSpeciID species"Arabidopsis thaliana"' % sys.argv[0]
        sys.exit(0)
    #------------------------------------------
    #------------------------------------------
    taxon_id_name = taxonID2Name()
    giTaxon = gi2taxon(sys.argv[2], taxon_id_name)
    giToSpeciID = giToID(sys.argv[3])
    species = sys.argv[4]
    #------------------------------------------
    mDict = {}
    group = 0
    for line in open(sys.argv[1]):
        if line[0] == '=':
            num = int(line.split()[1])
            if num:
                search = line[1:].rsplit('.',1)[0]
                group = 1
                mDict[search] = set()
            else:
                group = 0
        elif group and line[0] == '>':
            hit = line[1:-1]
            if hit.isdigit():
                hitTaxon = giTaxon[hit]
                ##hit slef species, detect if hit locus is search
                if hitTaxon == species: 
                    giagi = giToSpeciID[hit]
                    if search.find(giagi) != 0:
                        mDict[search].add(hitTaxon)
                ##hit no taxon, detect if it is species. If yes,
                ##detect if hit locus is search.
                elif hitTaxon == 'NoTaxon':
                    if hit in giToSpeciID:
                        giagi = giToSpeciID[hit]
                        if search.find(giagi) != 0:
                            mDict[search].add(species)
            #not all digits hit            
            else:
                if search.find(hit) != 0:
                    mDict[search].add(species)
        #----END one line----------------------------
    #--------END reading--------------------
    symbolDict = {0:'^', 1:'*'}
    for key, valueL in mDict.items():
        lenvalueL = len(valueL)
        removeone = 1 if species in valueL else 0
        if lenvalueL:
            print '>%s %d%s' % (key, lenvalueL, symbolDict[removeone])
            print '#'.join(valueL)
#------------END main-----------------


if __name__ == '__main__':
    main()

