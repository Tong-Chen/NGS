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
    parser.add_option("-o", "--output-folder", dest="output", default="scheme", 
        metavar="FILEIN", help="LOWER-CASE o. The directory to save generated output files. Default <scheme>.")
    parser.add_option("-O", "--override", dest="override",
        default=False, action="store_true", 
        help="Override last initialization. Upper-case O. Default no-override.")
    parser.add_option("-n", "--name-bookdown", dest="name",
        metavar="NAME", help="Name of bookdown file.")
    parser.add_option("-t", "--title", dest="title",
        metavar="TITLE", help="Title for document.")
    parser.add_option("-a", "--author-list", dest="author",
        metavar="AUTHOR", help="Author or affilication or other information. \
Multiple items separated by <;> can be given here. In format like \
客户单位：上帝;服务单位：易汉博基因科技（北京）有限公司;联系方式：vip@ehbio.com")
    parser.add_option("-d", "--document-class", dest="doc_class", default="article", 
        help="Document class, article,  book, report. Default [article]")
    parser.add_option("-l", "--latex-template", dest="latex_temp", 
        default="/home/ct/resource/ehbio.latex", 
        help="A file to set latex output template. [Optional, default /MPATHB/self/pandoc/ehbio.latex]")
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
def generateStyle_css(output, check=False):
    if check and os.path.exists(output+"/style.css"):
        print >>sys.stderr, "Looks like you have initiated the report. Pease make sure you do want to override."
        sys.exit(1)
    fh = open(output+"/style.css", 'w')
    print >>fh, '''
p.caption {
  color: #777;
  margin-top: 10px;
  text-align: justify;
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

body {
  text-align: justify;
}

.book .book-body .page-wrapper .page-inner section.normal table {
  display: block;
  overflow: auto;
  width: 100%;
}

.book .book-body .page-wrapper .page-inner section.normal caption {
  text-align: justify;
  width:100%;
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


def generateProj(output, name):
    file = name+'.Rproj'
    fh = open(output+'/'+file,'w')
    print >>fh,'''
Version: 1.0

RestoreWorkspace: Default
SaveWorkspace: Default
AlwaysSaveHistory: Default

EnableCodeIndexing: Yes
UseSpacesForTab: Yes
NumSpacesForTab: 2
Encoding: UTF-8

RnwWeave: Sweave
LaTeX: pdfLaTeX

BuildType: Website
'''
    fh.close()
#------------------------------------------
def formatAuthor(author, sep=';'):
    '''
    Transfer 'company1 - author1;company2 - author2; contact-infot' into

    - "company1 - author1"
    - "company2 - author2"
    - "contact-infot"
    '''
    tmpL = ["- \""+tmp.strip()+"\"" for tmp in author.split(sep)]
    return '\n'.join(tmpL)
#------------------------

def generateTestBib(output, bib):
    bibL = [i.strip() for i in bib.split(',')]
    for file in bibL:
        if not os.path.exists(file):
            fh = open(output+'/'+file, 'w')
            print >>fh, '''
@article{chen_m6a_2015,
	title = {m6A {RNA} {Methylation} {Is} {Regulated} by {MicroRNAs} and {Promotes} {Reprogramming} to {Pluripotency}},
	volume = {16},
	issn = {1934-5909, 1875-9777},
	url = {http://www.cell.com/cell-stem-cell/abstract/S1934-5909(15)00017-X},
	doi = {10.1016/j.stem.2015.01.016},
	language = {English},
	number = {3},
	urldate = {2016-12-08},
	journal = {Cell Stem Cell},
	author = {Chen, Tong and Hao, Ya-Juan and Zhang, Ying and Li, Miao-Miao and Wang, Meng and Han, Weifang and Wu, Yongsheng and Lv, Ying and Hao, Jie and Wang, Libin and Li, Ang and Yang, Ying and Jin, Kang-Xuan and Zhao, Xu and Li, Yuhuan and Ping, Xiao-Li and Lai, Wei-Yi and Wu, Li-Gang and Jiang, Guibin and Wang, Hai-Lin and Sang, Lisi and Wang, Xiu-Jie and Yang, Yun-Gui and Zhou, Qi},
	month = mar,
	year = {2015},
	pmid = {25683224},
	pages = {289--301},
}
'''
            fh.close()
        else:
            os.system("/bin/cp -f "+file+' '+output+'/')
#-------------------------------

def generateBib(output, bib, bib_style):
    '''
    if bib is given,  generate

    bibliography: [1.bib, 2.bib]
    biblio-style: bib_style
    link-citations: yes
    '''
    if bib:
        generateTestBib(output, bib)
        bibL = ['bibliography: ['+bib+']', 'biblio-style: '+bib_style, 'link-citations: yes']
        return '\n'+'\n'.join(bibL)
    else:
        return ''
#-----------------------------

def generateRproj(output, name):
    rproj = open(output+'/'+name+'.Rproj', 'w')
    print >>rproj, '''Version: 1.0

RestoreWorkspace: Default
SaveWorkspace: Yes
AlwaysSaveHistory: Default

EnableCodeIndexing: Yes
UseSpacesForTab: Yes
NumSpacesForTab: 2
Encoding: UTF-8

RnwWeave: knitr
LaTeX: pdfLaTeX

AutoAppendNewline: Yes
StripTrailingWhitespace: Yes

BuildType: Website
'''
    rproj.close()

#--------------------------------------------

def generateIndex(output, title, authorL, doc_class, bibOpt):
    #print >>sys.stderr, title
    #print >>sys.stderr, authorL
    #print >>sys.stderr, doc_class
    #print >>sys.stderr, bibOpt
    index = open(output+"/index.Rmd", 'w')
    print >>index, '''--- 
title: "{}"
author: 
{}
date: "`r Sys.Date()`"
documentclass: {}{}
site: bookdown::bookdown_site
---

```{{r setup, include=FALSE}}
library(knitr)
output <- opts_knit$get("rmarkdown.pandoc.to")
output <- "html"
html = FALSE
latex = FALSE
opts_chunk$set(echo = FALSE, out.width="100%", fig.align="center", fig.show="hold", warning=FALSE, message=FALSE)
if (output=="html") {{
	html = TRUE
}}
if (output=="latex") {{
	opts_chunk$set(out.width="95%", out.height='0.7\\\\textheight', out.extra='keepaspectratio', fig.pos='H')
	latex = TRUE
}}
knitr::opts_chunk$set(cache=TRUE, autodep=TRUE)
mtime <- function(files){{
  lapply(Sys.glob(files), function(x) file.info(x)$mtime)
}}
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

def generateBookdown(output, name):
    '''
    '''
    bookdown = open(output+"/_bookdown.yml", 'w')
    print >>bookdown, '''book_filename: "{}"
output_dir: "{}"
delete_merged_file: true
language:
  ui:
    chapter_name: ""'''.format(name, name)
    
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
def generatePics(output):
    picD = "/home/ct/resource/ehbio"
    os.system("/bin/cp -rf "+picD+ ' '+output+'/')

def generateMock(output):
    pass
#-----------------------------------------
def generateOtherFiles(output):
    cmd = "echo -e '# 课题目的 {{#project_aim}}\n\n' >{}/01-aim.Rmd".format(output)
    os.system(cmd)
    
    design = "{}/02-design.Rmd".format(output)
    design_fh = open(design, 'w')
    print >>design_fh, '''# 课题设计 {#project_design}

课题采用双因素正交实验设计 (Table \@ref(tab:design))，每种因素组合有10个样品。分别检测其基因和非编码RNA表达谱的变化，以探索肥胖对心梗发生的调控作用。

Table: (\#tab:design) 课题设计.

```{r design, results="asis"}
design <- ";未发心梗;心梗
体重正常;10;10
肥胖小鼠;10;10"

design_mat <- read.table(text=design, sep=";", header=T, row.names=1)
knitr::kable(design_mat, format="markdown")
```
'''
    
    design_fh.close()

    cmd = "echo -e '# 课题方案 {{#project_scheme}}\n\n~~Finished content~~' >{}/03-scheme.Rmd".format(output)
    os.system(cmd)

    cmd = "echo -e '# 结果概览 {{#project_result_summary}}\n\n' >{}/04-result.Rmd".format(output)
    os.system(cmd)
    
    cmd = "echo -e '# 课题讨论 {{#project_discussion}}\n\n' >{}/06-discussion.Rmd".format(output)
    os.system(cmd)
    
    company = open(output+'/09-company.Rmd', 'w')
    print >>company, '''# 公司简介 {#company_intro}

```{r company-html, eval=html, out.width="99%"}
knitr::include_graphics(c("ehbio/company_1.png", "ehbio/company_2.png"))
```

```{r company-pdf, eval=latex, out.width="99%", out.height='0.99\\\\textheight', out.extra='keepaspectratio'}
knitr::include_graphics(c("ehbio/company_1.png", "ehbio/company_2.png"))
```
'''
    company.close()

    ref = open(output+'/10-references.Rmd', 'w')
    print >>ref, '''`r if (knitr:::is_html_output()) '# References {-}'`
'''
    ref.close()

#-------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    output = options.output
    name = options.name
    title = options.title
    author = options.author
    authorL = formatAuthor(author)
    
    os.system("mkdir -p "+output)
    generatePics(output)
    doc_class = options.doc_class
    latex_temp = options.latex_temp
    bib = options.bib
    bib_style = options.bib_style
    bibOpt = generateBib(output, bib, bib_style)
    #print >>sys.stderr, bibOpt
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------

    generateStyle_css(output)

    generateIndex(output, title, authorL, doc_class, bibOpt)

    generateBookdown(output, name)

    generateOutputYML(output, latex_temp)

    generateOtherFiles(output)
    generateProj(output, name)
#-------------------------------------------------

if __name__ == '__main__':
    main()


