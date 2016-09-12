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

def readInterpro(interproFile, interproDict, locusL = []):
    '''
********************************************************
This read interproscan result to a dict.
********************************************************

FileFormat:

    AT2G01080.1	62CB99082D3204BD	231	HMMPfam	PF03168	LEA_2	108	208	5.699999999999995E-12	T	21-Sep-2010	IPR004864	Late embryogenesis abundant protein, group 2	

ResultFormat:

    interproDict = {
        locus:{{(start,end):PF04004 Leo1 start-end IPR004864 ..},}
    }
    
    '''
    limit = 0
    if len(locusL) > 0:
        limit = 1
    for line in open(interproFile):
        line = line.strip()
        aList = line.split('\t')
        locus = aList[0]
        if limit and (locus not in locusL):
            continue
        key = (int(aList[6]), int(aList[7]))
        pfamdesp = ' '.join(aList[4:6])
        interprodesp = ' '.join(aList[11:])
        pos = '-'.join((aList[6], aList[7]))
        desp = ' '.join((pfamdesp, interprodesp, pos))
        if locus in interproDict:
            if key not in interproDict[locus]:
                interproDict[locus][key] = desp
            else:
                print >>sys.stderr, 'Wrong in interproanno'
                sys.exit(0)
        else:
            interproDict[locus] = {}
            interproDict[locus][key] = desp
    #--------End for-------------------------------------------
#--------------------------------------------------------------

def readseq(seqfile, seqdict):
    '''
********************************************************
This read protein sequence file to a dict.
********************************************************
FileFormat:

    AT2G01080
    246
    ATHCGTAGHGYYYGYIGIGUII

ResultFormat:

    seqdict = {locus: seq,}
    '''
    i = 1
    for line in open(seqfile):
        if i == 1:
            locus = line.rstrip()
        elif i == 3:
            if locus in seqdict:
                print >>sys.stderr, "Duplicate sequence in seqfile"
                sys.exit(0)
            seqdict[locus] = line.rstrip()
            i = 1
            continue
        i += 1
    #---------------------------------------
#-------------------------------------------
def readSeq(seqfile, seqdict):
    readseq(seqfile, seqdict)

def readRep(repfile, repdict):
    '''
********************************************************
This reads repetition file.
********************************************************
FileFormat(No matter With or without the last #):
    Exact
    
    >AT1G62760.1
    SSLSPSS:25#SSLSPSS:51#SSLSPSS:100#SSLSPSS:128#
    SLSPSSPPP:52#SLSPSSPPP:62#SLSPSSPPP:101#SLSPSSPPP:111#
    
    Fuzzy

    >AT1G31870.1
    PSKPEPR:  32:   0#PSPEPNR: 159:   2#PSPEPAR: 267:   1#
    NISPPRR: 115:   0#DMSPPRR: 146:   2#DLSPPRR: 202:   1#DLSPPRR:
    231:   1#

    Also the cds sequence is OK
ResultFormat:
    repdict = 
    {
    AT1G62760.1: 
        [
            {(25, 31): 'SSLSPSS', (51,57): 'SSLSPSS'},
            {(52, 60): 'SLSPSSPPP',},    
        ] 
    }
'''
    locus = ''
    for line in open(repfile):
        if line[0] == '>':
            locus = line[1:].strip()
        else:
            if locus == '':
                print >>sys.stderr, 'Wrong format in repfile'
                sys.exit(1)
            if locus not in repdict:
                repdict[locus] = []
            line = line.rstrip()
            line = line.rstrip('#')
            lineL = line.split('#')
            tmpDict = {}
            for seqpos in lineL:
                ##################################20110712
                #Rewrite to process fuzzy matched repetitions
                ##################################20110712
                seqposL = seqpos.split(':')
                seq = seqposL[0]
                start = int(seqposL[1])
                #seq, start = seqpos.split(':')
                #start = int(start)
                end = start + len(seq) - 1
                tmpDict[(start, end)] = seq
            #-------End of for---------------------
            repdict[locus].append(tmpDict)
        #-----------End of else--------------------
    #---------------End of for---------------------
#---------------End readRep-----------------------
def outputRep(aDict, file):
    locusL = aDict.keys()
    locusL.sort()
    fh = open(file, 'w')
    for key in locusL:
        itemDL = aDict[key]
        #patched at 20110908, if there is no rep,no output
        if not len(itemDL):
            continue
        print >>fh, '>%s' % key
        for itemD in itemDL:
            groupL = []
            itemDKeyL = itemD.keys()
            itemDKeyL.sort()
            for posset in itemDKeyL:
                groupL.append(':'.join((itemD[posset],str(posset[0]))))
            #-------------------------------------------
            print >>fh, '#'.join(groupL)
        #-----------------------------------------------------
    #---------------------------------------------------------
    fh.close()
#----------------------------------------

def readseqrep(seqrepfile, seqdict, repdict):
    '''
********************************************************
This read protein sequence file with one more line of 
repetition or other sequence to a dict.
********************************************************
FileFormat:
    AT2G01080
    246
    ATHCGTAGHGYYYGYIGIGUII
    ATHCGTAGHGY:1#GYIGI:14(People Number)

ResultFormat:
    seqdict = {AT2G01080:"ATHCGTAGHGYYYGYIGIGUII",}
    repdict = {AT2G01080: [{(start, end):seq,}],}
    '''

    i = 0
    for line in open(seqrepfile):
        i += 1
        if i == 1:
            locus = line[:-1]
        elif i == 3:
            seq = line[:-1]
            if locus in seqdict:
                print >>sys.stderr, "Duplicate sequence in seqrepfile"
                sys.exit(0)
            seqdict[locus] = seq
        elif i == 4:
            i = 0
            replist = line[:-1].split('#')
            #print replist
            if locus in repdict:
                print >>sys.stderr, "Duplicate rep in seqrepfile"
                sys.exit(0)
            repdict[locus] = []
            tmpDict = {}
            for rep in replist:
                seq, start = rep.split(':')
                start = int(start)
                end = start + len(seq) -1
                key = (start, end)
                tmpDict[key] = seq
            repdict[locus].append(tmpDict)
            #---------------------------------------------                    
        #-------------------------------------------------
    #------------------------------------------------------
#----------------------------------------------------------
def readAnno(annofile, annodict, head = 1, locusL = []):
    '''
****************************************************************
This read anno file get from Ensmbl using R-Script.
****************************************************************
FileFormat:
    ensembl_peptide_id      description     go_biological_process_id    name_1006       go_cellular_component_idgo_cellular_component_id    go_cellular_component__dm_name_1006     go_molecular_function_id        go_molecular_function__dm_name_1006     interpro        interpro_description
    AT3G18710.1     plant U-box 29.[Source:TAIR;Acc:AT3G18710]      GO:0010200      response to chitin      GO:0000151      ubiquitin ligase complex        GO:0004842      ubiquitin-protein ligase activity       IPR003613       U box domain
    
ResultFormat:
    
    annodict = {
        AT2G01080 : plant U-box 29.\\GO:0004842 ubiquitin;
            GO:***\\Interpro ***** 
    }
    '''
    limit = 0
    if len(locusL) > 0:
        limit = 1
    tmpDict = {}
    for line in open(annofile):
        if head:
            head -= 1
        else:
            lineL = line[:-1].split('\t')
            key = lineL[0]
            if limit and (key not in locusL):
                    continue
            if key not in tmpDict:
                desp = set()
                go = set()
                interpro = set()
                tmpDict[key] = [desp, go, interpro]
            if lineL[1] != 'NA':
                tmp203 = lineL[1].replace(';', '; ')
                tmpDict[key][0].add(tmp203)
            if lineL[2] != 'NA' and lineL[3] != 'NA' and \
                len(lineL[2]) and len(lineL[3]):
                tmpDict[key][1].add(lineL[2]+ ': ' +lineL[3])
            if lineL[4] != 'NA' and lineL[5] != 'NA' and \
                len(lineL[4]) and len(lineL[5]):
                tmpDict[key][1].add(lineL[4]+ ': ' +lineL[5])
            if lineL[6] != 'NA' and lineL[7] != 'NA' and \
                len(lineL[6]) and len(lineL[7]):
                tmpDict[key][1].add(lineL[6]+ ': ' +lineL[7])
            if lineL[8] != 'NA' and lineL[9] != 'NA' and \
                len(lineL[8]) and len(lineL[9]):
                tmpDict[key][2].add(lineL[8]+ ": " +lineL[9])
        #----------End of else------------------------
    #--------------End of for-------------------------
    for key, value in tmpDict.items():
        annodict[key] = ''
        for ele in value:
            if len(ele):
                annodict[key] += '; '.join(ele) + '\\\\'
#-----------------------------------------------------

def readAnnoNew(annofile, annodict, head = 1, locusL = []):
    '''
****************************************************************
This read anno file get from Ensmbl using R-Script.
****************************************************************
FileFormat:
    ensembl_peptide_id      description     go_biological_process_id    name_1006       go_cellular_component_idgo_cellular_component_id    go_cellular_component__dm_name_1006     go_molecular_function_id        go_molecular_function__dm_name_1006     interpro        interpro_description
    ensembl_peptide_id	description	go_id	name_1006	definition_1006	namespace_1003	interpro	interpro_description	pfam
    AT3G18710.1     plant U-box 29.[Source:TAIR;Acc:AT3G18710]      GO:0010200      response to chitin      GO:0000151      ubiquitin ligase complex        GO:0004842      ubiquitin-protein ligase activity       IPR003613       U box domain
    AT1G01010.1	NAC domain containing protein 1.[Source:TAIR;Acc:AT1G01010]	GO:0007275	multicellular organismal development	"The biological process whose specific outcome is the progression of a multicellular organism over time from an initial condition (e.g. a zygote or a young adult) to a later condition (e.g. a multicellular animal or an aged adult)." [GOC:dph, GOC:ems, GOC:isa_complete, GOC:tb]	biological_process	IPR003441	No apical meristem (NAM) protein	PF02365

ResultFormat:
    
    annodict = {
        AT2G01080 : plant U-box 29.\\GO:0004842 ubiquitin;
            GO:***\\Interpro ***** 
    }
    '''
    limit = 0
    if len(locusL) > 0:
        limit = 1
    tmpDict = {}
    '''
    tmpDict = {locus: [despset("plant U-box 29"), 
            godict(
                {biological_process:(GO:10092 name \\\\definition) }
                {molecular_function:(GO:10092 name \\\\definition) }
                {cellular_component:(GO:10092 name \\\\definition) }
            ), 
            interproset("interproID interproDesp")
            pfamset("pfamid")
            ]}
    '''
    for line in open(annofile):
        if head:
            head -= 1
        else:
            (key, desp_2, go_3, gona_4, gode_5, goac_6, inter_7,
                interDes_8, pfam_9) = line[:-1].split('\t')
            if limit and (key not in locusL):
                    continue
            if key not in tmpDict:
                desp = set()
                godict = dict()
                interpro = set()
                pfam = set()
                tmpDict[key] = [desp, godict, interpro, pfam]
            if desp_2 != 'NA':
                tmp203 = desp_2.replace(';', '; ')
                tmpDict[key][0].add(tmp203)
            if go_3 != 'NA' and gona_4 != 'NA' and len(go_3) and \
                len(gona_4) and gode_5 != 'NA' and goac_6 != 'NA' and \
                len(gode_5) and len(goac_6):
                goac_6 = goac_6.replace('_', r'\_')
                if goac_6 not in tmpDict[key][1]: 
                    (tmpDict[key][1])[goac_6] = set()
                (tmpDict[key][1])[goac_6].add(go_3+" "+gona_4+ 
                    "\\\\"+gode_5+'\\\\')
                
            if inter_7 != 'NA' and interDes_8 != 'NA' and \
                    len(inter_7) and len(interDes_8): 
                tmpDict[key][2].add(inter_7 + "\t" + interDes_8)
            if pfam_9 != 'NA' and len(pfam_9):
                tmpDict[key][3].add(pfam_9)
        #----------End of else------------------------
    #--------------End of for-------------------------
    for key, value in tmpDict.items():
        if len(value[0]):
            annodict[key] = ' '.join(value[0]) + '\\\\'
        if len(value[1]):
            keyL = value[1].keys()
            keyL.sort()
            for keygoac in keyL:
                annodict[key] += keygoac+ '\\\\'+ \
                    '\\\\'.join(value[1][keygoac]) + '\\\\'
        if len(value[2]):
            annodict[key] += '; '.join(value[2]) + '\\\\'
        if len(value[3]):
            annodict[key] += '; '.join(value[3])
#-----------------------------------------------------
def readIsoformId(file):
    '''
    This read and find all the isoformers of one proteins.
    '''
    isoformDict = {}
    for line in open(file):
        line = line.strip()
        dot = line.find('.')
        assert dot != -1
        pre = line[:dot]
        after = line[dot:]
        if pre not in isoformDict:
            isoformDict[pre] = []
        isoformDict[pre].append(after)
    #----------------------------------
    return isoformDict

def readFasta(file, locusL=[]):
    '''
    >locus |annotation
    or
    >locus|annotation
    asdsdlsdkflsdsdsldlsdkfsdsd
    asdsdlsdkflsdsdsldlsdkfsdsd
    TAIR_CDS, PEP fills first format.
    SATIVA_CDSi, _PEP fills second 
    '''
    useL = 0
    if len(locusL) > 0:
        useL = 1
    #------------------------------------------
    aDict = {}
    for line in open(file):
        if line[0] == '>':
            saveThis = 1
            #changed 20110824
            #locus = (line[1:].split(' ')[0]).strip()
            locus = (line[1:].split('|')[0]).strip()
            if useL and locus not in locusL:
                saveThis = 0
            else:
                aDict[locus] = ''
        elif saveThis:
            #modified 20110824 to deal with seuqences with a '*' at
            #the end.
            #aDict[locus] += line.strip() #before
            line = line.strip()
            aDict[locus] += line.strip('*')
            #-------------------------------------------
    return aDict
#------------------readFasta----------------------------


def readStdFasta(file):
    '''
    >locus |annotation
    or
    >locus|annotation
    asdsdlsdkflsdsdsldlsdkfsdsd
    asdsdlsdkflsdsdsldlsdkfsdsd
    TAIR_CDS, PEP fills first format.
    SATIVA_CDSi, _PEP fills second 
    '''
    #------------------------------------------
    aDict = {}
    locus = ''
    for line in open(file):
        if line[0] == '>':
            if locus:
                aDict[locus] = ''.join(tmpL)
            locus = line[1:].strip()
            tmpL = []
        else:
            line = line.strip()
            tmpL.append(line.strip('*'))
            #-------------------------------------------
    #----------------------------------------------
    if locus:
        aDict[locus] = ''.join(tmpL)
    return aDict
#------------------readStdFasta----------------------------



def readStdPro(file, noduplex=1):
    '''
    >LOCUS|sth else
    FLMIVSPTAYHGNKDECWSRGPHELEUILEME
    TVALSERPROTHRALATYRHISGLNASNLYSA
    SPGLUCYSTRPARGSERARGGLY*
    '''
    aDict = {}
    for line in open(file):
        if line[0] == '>':
            locus = (line.split('|')[0])[1:]
            locus = locus.strip()
            if noduplex and locus in aDict:
                    print >>sys.stderr, 'Duplicate locus %s' % locus
            #---------------------------------------
            aDict[locus] = ''
        else:
            line = line.strip('\n*\r')
            aDict[locus] += line
    #---------------------------------------------------
    return aDict
#---------------------------------------------------------
def readColumns(file, head=1, sep='\t', first=0, second=1):
    '''
    This function reads in two columns and returns a dict,
    with <first>+1 column as key and <second>+1 column as value.
    <first> and <second> is the abnormal column index beginning at
    0. The columns are separated by '\t' defaultly. 
    '''
    aDict = {}
    for line in open(file):
        if head:
            head -= 1
            continue
        #-------------------------
        lineL = line.strip().split(sep)
        key, value = lineL[first], lineL[second]
        if key not in aDict:
            aDict[key] = set()
        aDict[key].add(value.upper())
    #-----------end reading------------
    return aDict
#----------------------------------------------------

def openCT(file):
    '''
    This aliases the default open in python.
    '''
    pass

#-----------------------------------------------------------


if __name__ == '__main__':
    print >>sys.stderr, '''
This is a package containing all the input and output functions which
can be called by other source code.
    '''

