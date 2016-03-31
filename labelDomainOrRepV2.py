#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
This is version 2.
1. Make the output of proteins in proper order if the locusfile are
notnot given. For the function dict.keys() is random, so sort
locusList. 
2. Modify the bug of loosing star after the first repetition when 
two repetitions are tandem and both have overlap with one domain.
Just like this:
D----------------------------------------
R1  -------------------=======================R2
3. Update label of repetitions with multiple stars.
4. Update explanation.
5. Alter annotation.
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
from optparse import OptionParser as OP
import sys
import ctIO
import os
#------------------common------------------------------------
numDict = {1:'one', 2:'two', 3:'three', 4:'four', 5:'five',
        6:'six', 7:'seven', 8:'eight', 9:'nine', 10:'ten', 11:'eleven'}
domainColDict = {1:'\'{Apricot}', 11:'\'{Mulberry}', 3:'\'{Salmon}', 
        4:'\'{SpringGreen}', 5:'\'{Gray}', 6:'\'{CornflowerBlue}',
        7:'\'{Tan}', 8:'\'{Cyan}', 9:'\'{Peach}',
        10:'\'{YellowGreen}', 2:'\'{SkyBlue}'}
repColDict = {1:'\'{Yellow}', 2:'\'{Melon}', 3:'\'{Green}',
        4:'\'{YellowOrange}', 5:'\'{GreenYellow}'}

texshade = r'''\seqtype{P}
\showsequencelogo{top}
\ttopspace{-\baselineskip}
\feature{top}{1}{1..300}{color:charge}{}
\showfeaturestylename{top}{charge}
\feature{bottom}{1}{1..300}{color:molweight[ColdHot]}{}
\showfeaturestylename{bottom}{weight}
\bbottomspace{-\baselineskip}
\feature{bbottom}{1}{1..300}{bar:hydrophobicity[Red,Gray10]}{}
\showfeaturestylename{bbottom}{hydrophob.}
\bargraphstretch{3}
\featurestylenamescolor{Red}
\featurestylenamesrm \featurestylenamesit
\hideconsensus
\end{texshade}
'''
#'''
#Special character:   #   $   %  ^       &   _   {    }   ~
#Transferred meaning: \# \$  \%  \^{}    \&  \_ \{   \}  \~{}
#'''
#------------------common------------------------------------


def cmdpara(argv):
    if len(argv) == 1:
        cmd = 'python '+ argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    if len(argv) == 3:
        print >>sys.stderr, 'Less parameters'
        sys.exit(1)
    desc = "Label domain or rep in mother sequence\n\n\
*********Version 2*******************\n"
    usages = desc+"%prog [-i interpro] [-a anno] [-s seq] \
[-r rep] [-f seqrep] [-l locus]"
    parser = OP(usage=usages)
    parser.add_option("-i", "--in", dest="interpro", metavar="INTERPRO",
        help="The output of intproscan")
    parser.add_option("-a", "--an", dest="anno", metavar="ANNO",
        help="The annotation get from ensembl")
    parser.add_option('-s', "--sq", dest="seqfile", metavar="SEQFILE",
        help="The protein sequence file")
    parser.add_option('-r', "--rp", dest="repfile", metavar="REPFILE",
        help="The repetition file")
    parser.add_option('-f', "--sr", dest="seqrepfile", metavar="SEQREPFILE",
        help="Put sequence and rep together")
    parser.add_option('-l', '--lc', dest="locusfile", metavar="LOCUSFILE",
        help="The locus file")
    parser.add_option('-o', '--order', dest='sort', default=False,
        help="Sort the locus[%default]", action="store_true")
    (options, args) = parser.parse_args(argv[1:])
    if options.interpro == None or options.anno == None:
        print >>sys.stderr, "Wrong parameters"
        sys.exit(1)
    return (options, args)
#---------cmdpara-----------------------------------------


def latexHead():
    print r'''\documentclass[a4paper, 9pt]{article}
\author{\href{mailto:chentong\_biology@163.com}{chentong}}
\title{My Work}
\date{\today}
\usepackage[top=0.3in, bottom=0.3in, left=0.3in, right=0.3in]{geometry}
\usepackage[colorlinks, bookmarks=true, pdfstartview=FitH, 
pdftitle=Protein repetition, pdfauthor=ct586]{hyperref}
\usepackage[usenames, dvipsnames]{color}
\usepackage{dnaseq}
\usepackage{texshade}
\newcommand{\anno}[1]{{\textcolor{Purple}{#1}}}
\newcommand{\repone}[1]{{\colorbox{Yellow}{#1}}}
\newcommand{\reptwo}[1]{{\colorbox{Melon}{#1}}}
\newcommand{\repthree}[1]{{\colorbox{Green}{#1}}}
\newcommand{\repfour}[1]{{\colorbox{YellowOrange}{#1}}}
\newcommand{\repfive}[1]{{\colorbox{GreenYellow}{#1}}}
\newcommand{\domainone}[1]{{\colorbox{Apricot}{#1}}}
\newcommand{\domaintwo}[1]{{\colorbox{SkyBlue}{#1}}}
\newcommand{\domainthree}[1]{{\colorbox{Salmon}{#1}}}
\newcommand{\domainfour}[1]{{\colorbox{SpringGreen}{#1}}}
\newcommand{\domainfive}[1]{{\colorbox{Gray}{#1}}}
\newcommand{\domainsix}[1]{{\colorbox{CornflowerBlue}{#1}}}
\newcommand{\domainseven}[1]{{\colorbox{Tan}{#1}}}
\newcommand{\domaineight}[1]{{\colorbox{Cyan}{#1}}}
\newcommand{\domainnine}[1]{{\colorbox{Peach}{#1}}}
\newcommand{\domainten}[1]{{\colorbox{YellowGreen}{#1}}}
\newcommand{\domaineleven}[1]{{\colorbox{Mulberry}{#1}}}
%\newcommand{\tair}[1]{\href{http://www.arabidopsis.org/servlets/Search?type=general\&search\_action=detail\&method=1\&show\_obsolete=F\&name=#1\&sub\_type=gene\&SEARCH\_EXACT=4\&SEARCH\_CONTAINS=1}{#1}}
\newcommand{\tair}[1]{\href{http://www.arabidopsis.org/servlets/TairObject?name=#1\&type=locus}{#1}}
\begin{document}
\maketitle
\tableofcontents
\clearpage
'''

def latexTail():
    print r'''
\end{document}    
'''

def latexExplain():

    print r'''\section{Intro}
Here is the description.

The functional domains are predicted by interproscan using HMM-Pfam.
This method can get most of the domains. You can refer to the link to
Arabidopsis web-site to see if there are more domains.

    '''
    for i in range(1, 6):
        print ''.join((r'\rep', numDict[i], '{Rep', str(i), '*'*i, '}'))
        print
    for i in range(1,12):
        print ''.join((r'\domain', numDict[i], '{Domain', str(i), '}'))
        print
    print r'\clearpage'
#-------------------------------------------------------------------



def getnewDictDomain(interSonDict, domainDespList, newDict):
    '''
    interSonDict:
        {(25,31): 'sdfsdfsdf',}
    '''
    domainPosL = interSonDict.keys()
    domainPosL.sort()
    domainNo = 1
    for domainPos in domainPosL:
        if domainNo > 11:
            domainNo = 11
        #if domainPos[0] not in newDict:
        desp = interSonDict[domainPos]
        desp = desp.replace('_', r'\_')
        desp = desp.replace('&', r'\&')
        newDict[domainPos[0]] = ''.join((r'\domain',\
            numDict[domainNo], '{'))
        lenDesp = len(desp)
        time = lenDesp / 90 + 1
        for eachT in range(time):
            end = (eachT+1) * 90 if (eachT+1)*90 < lenDesp \
                else lenDesp
            eachDesp = desp[eachT*90:end]
            domainDespList.append(''.join((newDict[domainPos[0]],\
                eachDesp, '}')))
        newDict[domainPos[0]] = domainColDict[domainNo]
        #else:
        #    print >>sys.stderr, 'Duplicate DomainPos'
        #    sys.exit(1)
        #if domainPos[1] not in newDict:
        newDict[domainPos[1]] = '\'{White}'
        #else:
        #   print >>sys.stderr, 'Duplicate DomainPos'
        #    sys.exit(1)
        domainNo += 1
#-----------------------------------------------------------

def getnewDictRepDesp(id, num, repDespList):
    '''
    repdictSonL
        [
            {(25, 31): 'SSLSPSS', (51,57): 'SSLSPSS'},
            {(52, 60): 'SLSPSSPPP',},    
        ] 
    '''
    for no in range(num):
        filename = id + str(no) + '.aln'  
        string = r'\begin{texshade}{tcoffee/' + filename + '}\n' + texshade
        repDespList.append(string)
        
#--------------------------------------------------------

def getnewDictRep(repdictSonL, domainPosL, repDespList, newDict):
    '''
    repdictSonL
        [
            {(25, 31): 'SSLSPSS', (51,57): 'SSLSPSS'},
            {(52, 60): 'SLSPSSPPP',},    
        ] 
    domainPosL
        [(25,31),()]
    '''
    noSonD = 0
    noCmp = 0
    if len(domainPosL) == 0:
        noCmp = 1
    #print >>sys.stderr,'begin', newDict, 
    for repdictSonD in repdictSonL:
        noSonD += 1
        if noSonD > 5:
            noSonD = 5
        repdictSonDK = repdictSonD.keys()
        repdictSonDK.sort()

        for repPosT in repdictSonDK:
            #---------Here for desp----------------------------------
            #repDespList.append(''.join((r'\rep', numDict[noSonD], '{',\
            #    repdictSonD[repPosT], str(repPosT), '}')))
            #----------------------------------------------------------
            #----------Here for label----------------------
            #print >>sys.stderr, 'New', newDict
            if noCmp:
                newDict[repPosT[0]] = repColDict[noSonD]
                newDict[repPosT[1]] = '*' * noSonD + '\'{White}'
                continue
            overlap = 0
            for domainPosT in domainPosL:
                #-------Consider pre number----------------------------
                #print >>sys.stderr, 'Has domain', newDict
                #print >>sys.stderr, domainPosT
                #print >>sys.stderr, repPosT
                if repPosT[0] == domainPosT[0]:
                    if repPosT[1] < domainPosT[1]:
                        #D--------------------
                        #R-----------
                        old = newDict[domainPosT[0]]
                        if repPosT[0] in newDict:
                            newDict[repPosT[0]] += repColDict[noSonD]
                        else:
                            newDict[repPosT[0]] = repColDict[noSonD]
                        newDict[repPosT[1]+1] = old + '*' * noSonD
                        overlap = 1
                    elif repPosT[1] > domainPosT[1]:
                        #D-------
                        #R-----------
                        #del newDict[domainPosT[1]]
                        newDict[domainPosT[1]+1] = repColDict[noSonD]
                        newDict[repPosT[1]] = '\'{White}' + '*' * noSonD 
                        overlap = 1
                    else:
                        #D-------
                        #R-------
                        newDict[repPosT[1]] = \
                            newDict[repPosT[1]] + '*' * noSonD
                        overlap = 1
                #---------left equal----------------------------
                elif domainPosT[0] < repPosT[0] < domainPosT[1]:
                    newDict[repPosT[0]] = repColDict[noSonD]
                    if repPosT[1] < domainPosT[1]:
                        #D-----------------------
                        #R       -----------
                        old = newDict[domainPosT[0]]
                        newDict[repPosT[1]+1] = old + '*' * noSonD
                        overlap = 1
                    elif repPosT[1] == domainPosT[1]:
                        #D-----------------------
                        #R       ----------------
                        newDict[repPosT[1]] = \
                            newDict[domainPosT[1]]+'*' * noSonD
                        overlap = 1
                    else:
                        #D-----------------------
                        #R       ----------------------
                        #del newDict[domainPosT[1]]
                        newDict[domainPosT[1]+1] = repColDict[noSonD]
                        newDict[repPosT[1]] = '\'{White}'+'*' * noSonD 
                        overlap = 1
                #-------left in-------------------------------
                elif repPosT[0] < domainPosT[0]:
                    newDict[repPosT[0]] = repColDict[noSonD]
                    if domainPosT[0] < repPosT[1] < domainPosT[1]:
                        #D      ------------------
                        #R-------------
                        old = newDict[domainPosT[0]]
                        newDict[repPosT[1]+1] = old + '*' * noSonD
                        #del newDict[domainPosT[0]]
                        newDict[domainPosT[0]+1] = repColDict[noSonD]
                        overlap = 1
                    elif repPosT[1] == domainPosT[1]:
                        #D      -------
                        #R-------------
                        newDict[repPosT[1]] = \
                            newDict[domainPosT[1]] + '*' * noSonD
                        overlap = 1
                    elif domainPosT[1] < repPosT[1]:
                        #D      -------
                        #R------------------
                        #del newDict[domainPosT[1]]
                        newDict[domainPosT[1]+1] = repColDict[noSonD]
                        newDict[repPosT[1]] = '\'{White}' + '*' * noSonD 
                        overlap = 1
                #others
                #D------          ----------
                #R  ------------------
                #if overlap:
                #    break
                #----left out--------------------------------------
                #-------Consider pre number----------------------------
            #---Compare domain rep pos-------------------------------------
            if overlap == 0:
                newDict[repPosT[0]] = repColDict[noSonD]
                newDict[repPosT[1]] = '\'{White}' + '*' * noSonD 
        #----------repPosT-----------------------------------------
    #------------repdictSonD---------------------------------------
    #print >>sys.stderr, 'end', newDict
#----------------getnewDictRep-------------------------

def modifySeq(seq, newDict):
    newDictKeyL = newDict.keys()
    newDictKeyL.sort()
    shift = 0
    for newDictK in newDictKeyL:
        newDictV = newDict[newDictK] 
        if newDictV.find('White') != -1:
            seq.insert(newDictK+shift, newDictV)  
        else:
            seq.insert(newDictK+shift-1, newDictV)  
        shift += 1
#--------------------------------------------------------------------

def main():
    (options, args) = cmdpara(sys.argv)
    if options.sort:
        print 'sort'
        sys.exit(1)
    else:
        print 'no sort'
        sys.exit(1)
    print >>sys.stderr, "*******Print the result to screen.*******"

    #-------------------macro-------------------------------------
    isRep = 0
    isLoc = 0
    #-------------------------------------------------------------
    if options.seqfile != None:
        seqdict = {}
        ctIO.readseq(options.seqfile, seqdict)
    if options.seqrepfile != None:
        seqdict = {}
        repdict = {}
        ctIO.readseqrep(options.seqrepfile, seqdict, repdict)
        isRep = 1
    if options.repfile != None:
        repdict = {}
        ctIO.readRep(options.repfile, repdict)
        isRep = 1
    if options.locusfile != None:
        locusList = [line.strip() for line in open(options.locusfile)]
        isLoc = 1
    if not isLoc:
        locusList = repdict.keys() if isRep else seqdict.keys()
        locusList.sort()
    if isRep or isLoc:
        annodict = {}
        ctIO.readAnno(options.anno, annodict, 1, locusList)
        interproDict = {}
        ctIO.readInterpro(options.interpro, interproDict, locusList)

    #print locusList
    #print repdict.keys()
    #print repdict
    #print isRep
    #sys.exit(1)
    #----------------------------------------------------------------
    latexHead()
    latexExplain()
    #---------------------------------------------------------------
    for id in locusList:
        if id not in seqdict:
            print >>sys.stderr, "Unknown locus %s" % id
        else:
            hasInterpro = 0
            seq = list(seqdict[id])
            #-------------get newDict--------------------------------
            newDict = {}
            if id in interproDict:
                hasInterpro = 1
                domainDespList = []
                domainPosL = interproDict[id].keys()
                getnewDictDomain(interproDict[id], domainDespList,
                    newDict)
            #------------------------------------------------
            if isRep:
                repDespList = []
                repdictSonL = repdict[id]
                num = len(repdictSonL)
                if not hasInterpro:
                    domainPosL = []
                getnewDictRep(repdictSonL, domainPosL, repDespList, newDict)
                getnewDictRepDesp(id, num, repDespList)
            #----------------get newDict---------------------------           
            modifySeq(seq, newDict)
            #---------------------------------------------------------
            shortAnno = ''
            annos = annodict[id].replace('_', r'\_')
            annos = annos.replace('%', r'\%')
            annos = annos.replace('~', r'\~')
            annos = annos.replace('&', r'\&')
            firstBr = annos.find('[')
            if firstBr != -1:
                shortAnno = annos[:firstBr]
            print ''.join((r'\section{', id, ' ', shortAnno, '}' ))
            annos = r'\tair{' + id[:-2] + '} ' + annos
            print r'\anno{', annos, '}'

            print r'''
\noindent\begin{minipage}{\textwidth}
\noindent\rule{\textwidth}{2pt}
\DNA!'''
            #--without annotation
            seq = ''.join(seq)
            print seq

            print r'''!
\end{minipage}            
'''
            #------------------------------------------------
            print
            print '.' * 100
            print
            #---------------Rep------------------------------
            if isRep:
                for repSeq in repDespList:
                    print repSeq
                    print
            #-------------Domain desp----------------------
            if hasInterpro:
                print '.' * 100
                print
                for domainDesp in domainDespList:
                    print domainDesp
                    print
            print r'\clearpage'
            print
        #----------------End of else ---one locus-------------
    #----------------END of for ---all locus------------------
    latexTail()
    
    

if __name__ == '__main__':
    main()

