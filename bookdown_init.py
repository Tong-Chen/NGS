#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2016, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to generate the minimal requirement for bookdown usage.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
reload(sys)
sys.setdefaultencoding('utf8')

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
    parser.add_option("-o", "--output-folder", dest="output",
        metavar="FILEIN", help="The directory to save output files.")
    parser.add_option("-t", "--title", dest="title",
        metavar="TITLE", help="Title for document.")
    parser.add_option("-a", "--author-list", dest="author",
        metavar="AUTHOR", help="Author or affilication or other information. \
            Multiple items separated by <;> can be given here.")
    parser.add_option("-d", "--document-class", dest="doc_class", default="article", 
        help="Document class, article,  book, report. Default [article]")
    parser.add_option("-l", "--latex-template", dest="latex_temp", default="null", 
        help="A file to set latex output template. [Optional]")
    parser.add_option("-b", "--bibliography", dest="bib",
        help="Bibliography files for reference citation. Nomally these files can be got from EndNote or Zotero. Multiple files should be supplied as <1.bib, 2.bib>.")
    parser.add_option("-B", "--biblio-style", dest="bib_style", default="apalike", 
        help="Bibliography style for reference citation. Default <apalike>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.output != None, "A filename needed for -o"
    return (options, args)
#--------------------------------------------------------------------
def generateStyle_css(output):
    fh = open(output+"/style.css", 'w')
    print >>fh, '''
p.caption {
  color: #777;
  margin-top: 10px;
}
p code {
  white-space: inherit;
}
pre {
  word-break: normal;
  word-wrap: normal;
}
pre code {
  white-space: inherit;
}

.book .book-body .page-wrapper .page-inner section.normal table th {
  font-size: 70%;
}

.book .book-body .page-wrapper .page-inner section.normal table td {
  font-size: 70%;
}

.book .book-body .page-wrapper .page-inner section.normal ol {
  padding: 0 0 0 1em;
}
'''
#-----------------------------

def formatAuthor(author, sep=';'):
    '''
    Transfer 'company1 - author1;company2 - author2; contact-infot' into

    - "company1 - author1"
    - "company2 - author2"
    - "contact-infot"
    '''
    tmpL = ["- \""+tmp+"\"" for tmp in author.split(sep)]
    return '\n'.join(tmpL)
#------------------------

def generateBib(bib, bib_style):
    '''
    if bib is given,  generate

    bibliography: [1.bib, 2.bib]
    biblio-style: bib_style
    link-citations: yes
    '''
    if bib:
        bibL = ['bibliography: ['+bib+']', 'biblio-style: '+bib_style, 'link-citations: yes']
        return '\n'+'\n'.join(bibL)
    else:
        return ''
#-----------------------------

def generateIndex(output, title, authorL, doc_class, bibOpt):
    index = open(output+"/index.Rmd", 'w')
    print >>index, '''--- 
title: "{}"
author: 
{}
date: "`r Sys.Date()`"
documentclass: {}{}
---

```{{r setup, include=FALSE}}
library(knitr)
output <- opts_knit$get("rmarkdown.pandoc.to")
html = FALSE
latex = FALSE
opts_chunk$set(echo = FALSE, out.width="100%", fig.align="center", fig.show="hold")
if (output=="html") {{
	html = TRUE
}}
if (output=="latex") {{
	opts_chunk$set(out.width="95%", out.height='0.7\\textheight', out.extra='keepaspectratio', fig.pos='H')
	latex = FALSE
}}
#knitr::opts_chunk$set(cache=FALSE, autodep=TRUE)
set.seed(0304)
```

```{{asis, echo=html}}

# EHBIO Gene Technology {{-}}

```

```{{r cover, eval=html, out.width="99%"}}
knitr::include_graphics("ehbio/cover.png")

```
'''.format(title, authorL, doc_class, bibOpt)

    index.close()
#-------------------------------

def generateBookdown(output):
    '''
    '''
    bookdown = open(output+"/_bookdown.yml", 'w')
    print >>bookdown, '''book_filename: "{}"
output_dir: "{}"
language:
  ui:
    chapter_name: ""'''.format(output, output)
    
    bookdown.close()
#-----------------------------------

def generateOutputYML(output, latex_temp):
    output_yml = open(output+'/_output.yml', 'w')
    print >>output_yml, '''bookdown::pdf_book:
  template: {}
  latex_engine: xelatex
  citation_package: natbib
  keep_tex: yes
  toc_unnumbered: no
  toc_depth: 3
bookdown::epub_book: default
bookdown::gitbook:
  css: style.css
  config:
    download: [pdf]
    toc:
      before: |
          <li><a href="http://www.ehbio.com"><img src="http://www.ehbio.com/logos/ehbio_gitbook_logo.png" width="95%"></a></li>
      after: |
          <li><a href="mailto:ct@ehbio.com" target="blank">ct@ehbio.com</a></li>
    sharing: 
      github: no
      facebook: no
      twitter: no'''.format(latex_temp)
#------------------------------------------
def generateMock(output):
    pass
#-----------------------------------------
def generateOtherFiles(output):
    cmd = "echo '# 课题目的 {{#project_aim}}' >{}/01-aim.Rmd".format(output)
    os.system(cmd)

    cmd = "echo '# 课题设计 {{#project_design}}' >{}/02-design.Rmd".format(output)
    os.system(cmd)

    cmd = "echo '# 课题方案 {{#project_scheme}}' >{}/03-scheme.Rmd".format(output)
    os.system(cmd)

    cmd = "echo '# 课题结果 {{#project_result}}' >{}/04-result.Rmd".format(output)
    os.system(cmd)
    
    cmd = "echo '# 课题讨论 {{#project_discussion}}' >{}/04-discussion.Rmd".format(output)
    os.system(cmd)

    ref = open(output+'/10-references.Rmd', 'w')
    print >>ref, '''`r if (knitr:::is_html_output()) '# References {-}'`
'''
    ref.close()

#-------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    output = options.output
    title = options.title
    author = options.author
    authorL = formatAuthor(author)
    
    doc_class = options.doc_class
    latex_temp = options.latex_temp
    bib = options.bib
    bib_style = options.bib_style
    bibOpt = generateBib(bib, bib_style)
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    os.system("mkdir -p "+output)

    generateStyle_css(output)

    generateIndex(output, title, authorL, doc_class, bibOpt)

    generateBookdown(output)

    generateOutputYML(output, latex_temp)

    generateOtherFiles(output)
#-------------------------------------------------

if __name__ == '__main__':
    main()


