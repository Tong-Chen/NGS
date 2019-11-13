#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division,  with_statement

import sys
reload(sys)
sys.setdefaultencoding('utf8')

import urllib
from time import sleep

from mechanize import Browser
from bs4 import BeautifulSoup
br = Browser()
br.addheaders = [('User-agent',  'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1) Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]

import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

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
        #global desc
        #print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="Latin words of species or chinese names of species with one at a row.")
    parser.add_option("-t", "--type", dest="type",
        help="<latin> or <chinese> names")
    parser.add_option("-o", "--output", dest="output",
        help="Output file name")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def parseOneWord(species, type, output_fh):
    sleep(2)
    #print species
    if type == 'latin':
        type = '1'
    elif type == 'chinese':
        type = '2'
    #url_base = "http://www.cvh.org.cn"
    url_base = "http://www.cvh.ac.cn"
    url_search = url_base+"/search/"+urllib.quote(species.encode('utf8'))
    url = url_search+"?page=1&searchtype=1&n="+type

    nite_A = br.open(url)
    nite_A_html = BeautifulSoup(nite_A, 'html.parser')

    ## Get number

    # From em attribute

    em_words = nite_A_html.find(style="color:red")
    em_words_parent_div = em_words.parent
    em_words_parent_div_text = em_words_parent_div.text
    number_start = em_words_parent_div_text.find(u'共')+1
    number_end   = em_words_parent_div_text.find(u'号')
    number = int(em_words_parent_div_text[number_start:number_end])
    if number == 0:
        return 0
    pagenumber = int(number/15)+1+1

    
    parseOnepage(url_base, nite_A_html, 1, output_fh)
    ## Get link
    
    for page in range(2, pagenumber):
        url = url_search+"?page="+str(page)+"&searchtype=1&n="+type
        nite_A = br.open(url)
        nite_A_html = BeautifulSoup(nite_A, 'html.parser')
        parseOnepage(url_base, nite_A_html, page, output_fh)
        if page % 10 == 1:
            sleep(5)
        if page % 20 == 1:
            sleep(40)
        if page % 40 == 1:
            sleep(60)
        
#--------------------------
    
def parseOnepage(url_base, nite_A_html, page, output_fh):
    eachTDl = nite_A_html.find_all("td", attrs={"class":"td3"})
    eachTDl_len = len(eachTDl)
    i = 0

    for i in range(0, eachTDl_len, 7):
        eachTD_link = eachTDl[i].a.get('href')
        latin_name  = eachTDl[i+1].text
        chinese_name  = eachTDl[i+2].text
        sub_url = url_base + urllib.quote(eachTD_link)
        sub_nite_A = br.open(sub_url)
        sub_nite_A_html = BeautifulSoup(sub_nite_A, 'html.parser')
        sub_div = sub_nite_A_html.find_all("div", attrs={"id":"pe_qianc_2"})[0]

        if page == 1 and i == 0:
            tmpHeaderL = sub_div.find_all("span", attrs={"class":"sptitle2"})
            headerL = [span.string for span in tmpHeaderL]
            print >>output_fh, "Latin_name\tChinese\t{}".format("\t".join(headerL))
        #------------------------------------------
        #if i == 0:
        #    print "# {} {}".format(latin_name, chinese_name) 
        valueL = [latin_name, chinese_name]
        valueL.append(sub_div.find("div", id="o_spcollter").text)
        valueL.append(sub_div.find("div", id="o_spcoldate").text)
        valueL.append(sub_div.find("div", id="o_spplace").text)
        valueL.append(sub_div.find("div", id="o_spenviro").text)
        valueL.append(sub_div.find("div", id="o_spal").text)
        valueL.append(sub_div.find("div", id="o_sphabit").text)
        valueL.append(sub_div.find("div", id="o_sppreparations").text)
        valueL.append("")
        newValueL = []
        for empV in valueL:
            if not empV:
                empV='-'
            newValueL.append(empV)
        print >>output_fh, "\t".join(newValueL) 
    #-------------------------
#--------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    type = options.type
    output = options.output
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    output_fh = open(output, 'w')
    for line in fh:
        species = line.strip()
        parseOneWord(species, type, output_fh)
    output_fh.close()
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    #fh = open(output+'.finish', 'w')
    #print >>fh, "FINISHED"
    #fh.close()
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
    #if verbose:
    #    print >>sys.stderr,\
    #        "--Successful %s" % strftime(timeformat, localtime())
    #return 0

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'w')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()
    #return 0




