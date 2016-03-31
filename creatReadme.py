#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
This is a file to transfer README file to tex syntax.
You could operate it through several options.

Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__version__ = '0.1'
__revision__ = '0.1'
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
#from __future__ import division, with_statement
import sys

def header(texname, projectName):
    '''
    header(texname)

    this function is used to get a header file containing filename,
    author, time and other information.
    Then set the documentclass, the author information and the 
    projectName for the document.
    '''
    import time
    now = time.localtime(time.time())
    currentTime = time.strftime("%Y-%m-%d %H:%M:%S %a", now)
    print '%%        File: %s\n%%     Created: %s\n%%Last Changed: \
            %s\n%%' % (texname, currentTime, currentTime)
    print r'''
\documentclass[a4paper, 12pt]{article}
%define the title
\author{\href{mailto:chentong\_biology@163.com}{chentong}}'''
    print '\\title{%s}\n\\date{\\today}\n' % projectName


def usingListings():
    print r'''%this is used to display source code
\usepackage{listings}
\usepackage{xcolor}
%set parameter
\lstset{numbers=left, numberstyle=\tiny, keywordstyle=\color{blue!70},
commentstyle=\color{red!50!green!50!blue!50}, frame=shadowbox,
rulesepcolor=\color{red!20!green!20!blue!20}, escapeinside=``,
%chinese character should be in ``
xleftmargin=2em, xrightmargin=2em, aboveskip=1em}
%---------------------------------------listings----------------------
'''

#TODO we owe a titile here
def usingHyperlink():
    print r'''%---------------hyperlink------------------------------
\usepackage[colorlinks, bookmarks=true, pdfstartview=FitH, 
pdftitle=Protein repetition, pdfauthor=ct586]{hyperref}
%---------------hyperlink--------------------------------------------
    '''

def usingModule():
    usingListings()
    usingHyperlink()

def usingNewCommand():
    '''
    This is the using of self constructed command.
    '''
    print r'''
%--------------self constructed command------------------------------
\newcommand{\ctfoldlink}[1]{\href{http://210.75.224.29/ftp/project1/%
    #1}{#1}}
\newcommand{\descripfont}[1]{\small {#1}}
\newcommand{\fuzzyhref}[1]{\href{http://210.75.224.29/ftp/project1/%
    result/fuzzymatch/qualityEvaluate/#1}{#1}}
\newcommand{\exacthref}[1]{\href{http://210.75.224.29/ftp/project1/%
    result/exactmatch/qualityEvaluate/#1}{#1}}
\newcommand{\numcol}[1]{{\color{blue}{#1}}}
\newcommand{\report}[1]{\subsection{#1}\label{#1}\verbatimtabinput%
    [4]{document/#1}}
\newcommand{\findrepC}[1]{\subsection{#1}\label{#1}\lstinputlisting%
    [language=C]{../../cbin/findRep/#1}}
%--------------self constructed command------------------------------    
    '''

def beginDoc():
    '''
    From this, we begin the document.
    Using a directory to display the outline.
    You can choose whether you want a new page fot the main content.
    '''
    print r'''
\begin{document}
\maketitle
\tableofcontents
\clearpage
    '''

def parseDoc(inputfile):
    '''
    This is used to parse the input file to special format and add
    a tex label for each of them to make it can be read by latex.

    The file format it can discrimate is like the following.
    Ignore '(' and ')', you may want to substitue the content in 
    them to your words.
    The words after '#' is comment.
    All the website should be put in a {}, not allowed to be splited 
    in two lines.
    1.(Introduction)
    (This belongs to the project 'phosphoBlast'.)
        #the folder name: what this folder is used for(for folder name, 
        better no use '_', but whatever, now is ok)
        (phosphoblast): (this is a working folder.)
            #the file name: what this file is used for.(this can be a 
            subfolder, too)
            (abandonedIn3.75): this is the abandoned locus in version3.75
            which still in pepfile, total number is 48.
            # a little limit, this second line is not allowed to use colon.
    '''
    import re

    section = re.compile('[1-9]+\.(.+)')
    secDesc = re.compile('\S')
    listCon = '(\S+):\s*?(.*)'
    list0 = list1 = list2 = list3 = list4 = list5 = list6 =\
            list7 = list8 = list9 = ''
    listDesc0 = listDesc1 = listDesc2 = listDesc3 = listDesc4 = \
            listDesc5 = listDesc6 = listDesc7 = listDesc8 = \
            listDesc9 = ''
    listList = [list0, list1, list2, list3, list4, list5, list6,\
            list7, list8, list9]
    listDescList = [listDesc0, listDesc1, listDesc2, listDesc3,\
            listDesc4, listDesc5, listDesc6, listDesc7,\
            listDesc8, listDesc9] 
    for i in range(9):
        listCon = '\t' * i + listCon
        listList[i] = re.compile(listCon)
        listDescList[i] = re.compile('\t' * i + '(.+)')
    ##list = re.compile('\t(\S+):\s*?(.*)')
    ##listDesc = re.compile('\t(.+)')
    ##nlist = re.compile('\t\t(\S+):\s*?(.*)')
    ##nlistDesc = re.compile('\t\t(.+)')
    ##nnlist = re.compile('\t\t\t(\S+):\s*?(.*)')
    ##nnlistDesc = re.compile('\t\t\t(.+)')
    specialChar = re.compile('([_])')
    href = re.compile('({.+?})')
    #flags 
    listFlag = 0
    nlistFlag = 0
    nnlistFlag = 0

    for line in open(inputfile):
        line = specialChar.sub(r'\\\1', line)   
        line = href.sub(r'\href\1{here}', line)
        match = section.match(line)
        if match:
            if listFlag == 1:
                listFlag = 0
                print r'\end{enumerate}'
                print 
            raw = r'\section{'
            print '%s%s}' % (raw, match.group(1))
        else:
            match = secDesc.match(line)
            if match:
                print line,
            else:
                match = list0.match(line)
                if match:
                    if listFlag == 0:
                        listFlag = 1
                        print r'\begin{enumerate}'
                    if nlistFlag == 1:
                        nlistFlag = 0
                        print r'    \end{itemize}'
                    raw1 = r'  \item \ctfoldlink{'
                    print '%s%s}: %s' % (raw1, match.group(1), \
                            match.group(2))
                else:
                    match = list1.match(line)
                    if match:
                        if nlistFlag == 0:
                            nlistFlag = 1
                            print r'    \begin{itemize}'
                        if nnlistFlag == 1:
                            nnlistFlag = 0
                            print r'      \end{itemize}'
                        raw2 = r'      \item '
                        print '%s%s: %s' % (raw2,match.group(1),\
                                match.group(2))
                    else:
                        match = list3.match(line)
                        if match:
                            if nnlistFlag == 0:
                                nnlistFlag == 1
                                print r'      \begin{itemize}'
                            raw2 = r'        \item '
                            print '%s%s: %s'%(raw2,match.group(1),\
                                    match.group(2))
                        else:
                            sum = listFlag + nlistFlag + nnlistFlag
                            match = listDescList[sum-1].match(line)
                            ##if sum == 1:
                            ##    match = listDescList[sum-1].match(line)
                            ##elif sum == 2: 
                            ##   match = nlistDesc.match(line)
                            ##elif sum == 3:
                            ##   match = nnlistDesc.match(line)
                            if match:
                                space = '      ' * sum
                                print '%s%s' % (space, match.group(1))
    if listFlag == 1:
        print r'\end{enumerate}'

def endDoc():
    print
    print r'\end{document}'



def yieldTex(texname, projectName, inputfile):
    '''
    
    '''
    header(texname, projectName)
    usingModule()
    usingNewCommand()
    beginDoc()
    parseDoc(inputfile)
    
    endDoc()

def compileLatex(texname):
    '''
    This function is used to compile the created tex files.
    '''
    print 'Begin compiling'
    import os
    os.system("xelatex " + texname)
    print 'End compiling'

def main():
    if len(sys.argv) != 4:
        print 'Using python %s texname projectName inputfile' \
                % sys.argv[0]
        sys.exit(0)
    oldstdout = sys.stdout
    fh = open(sys.argv[1], 'w')
    sys.stdout = fh

    yieldTex(sys.argv[1], sys.argv[2], sys.argv[3])
    
    fh.close()
    sys.stdout = oldstdout
    print 'End'

    while True:
        print '''Please choose a command.
        c-->compile, q-->quit.
        If you add a new section, you may need to execute this more 
        than once.
        '''
        cmd = raw_input('>>>>')
        if cmd == 'c':
            compileLatex(sys.argv[1])
        elif cmd == 'q':
            break


if __name__ == '__main__':
    main()
