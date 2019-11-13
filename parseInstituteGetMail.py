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
    This is designed to parse institute pages and get tutors information.
    Current support IGDB, IOZ
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime, sleep
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import urllib
#from multiprocessing.dummy import Pool as ThreadPool

from mechanize import Browser
from bs4 import BeautifulSoup
br = Browser()
br.addheaders = [('User-agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1) Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]
br.set_handle_robots(False)

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
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="comma separated list of institutes. Currently support IGDB, IOZ, IB")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def parse_cksp_eol_cn(url):
    #http://cksp.eol.cn/tutor_search.php?do=q_search&keyword=%E7%94%9F%E7%89%A9&Submit=%E6%90%9C+%E7%B4%A2&key=research
    #http://cksp.eol.cn/tutor_search.php?page=1&do=q_search&keyword=%E7%94%9F%E7%89%A9&key=research
    # 35: url = "http://cksp.eol.cn/tutor_search.php?page="+str(i)+"&do=q_search&keyword=%E7%94%9F%E7%89%A9&key=research"
    # 2: url="http://cksp.eol.cn/tutor_search.php?page="+str(i)+"&do=q_search&keyword=%E5%8F%91%E8%82%B2&key=research"
    keywordL = ["肿瘤", "发育", "动物", "疾病", "生理", "遗传", "医学", "心脏", "基因", "神经", "临床", "免疫", "蛋白"]
    sub_urlS = set()
    for keyword in keywordL:
        keyword = urllib.quote(keyword.encode('utf8'))
        url = "http://cksp.eol.cn/tutor_search.php?do=q_search&keyword="+keyword+"&Submit=%E6%90%9C+%E7%B4%A2&key=research"
        #url = "http://cksp.eol.cn/tutor_search.php?page="+str(i)+"&do=q_search&keyword=%E5%8C%BB%E5%AD%A6&key=research"
        html = BeautifulSoup(br.open(url), 'html.parser')
        count = html.find('table', attrs={"class":"borBlue", "width":"666", "cellspacing":"5"})
        count = int(count.find_all("span")[-1].get_text())
        if not count:
            continue
        pages = int(count/21)+2
        for i in range(1, pages):
            if i > 1:
                url = "http://cksp.eol.cn/tutor_search.php?page="+str(i)+"&do=q_search&keyword="+keyword+"&key=research"
                html = BeautifulSoup(br.open(url), 'html.parser')
            #--------------------------------------------    
            html_table = html.find("table", attrs={'class':'tab_01'})
            html_table_tr = html_table.find_all("tr")
            for tr in html_table_tr[1:]:
                name = tr.find('td').get_text()
                sub_url = tr.find('td').find('a').get('href')
                if sub_url in sub_urlS:
                    continue
                sub_urlS.add(sub_url)
                html2 = BeautifulSoup(br.open(sub_url), 'html.parser')
                html2_table = html2.find("table", class_="tab_01")
                mail = html2_table.find_all('td')[1].get_text()
                print '\t'.join([name, mail, '-', sub_url])
                sleep(2)
            sleep(20)
        sleep(40)
#--------------------------------------------------

def parseEachTutorIGDB(url):
    html = BeautifulSoup(br.open(url), 'html.parser')
    tutor_name = html.find('h1')
    if tutor_name:
        tutor_name = tutor_name.get_text()
    else:
        tutor_name = '-'
    phone = html.find('div', id="dh")
    if phone:
        phone = phone.get_text()
    else:
        phone = '-'
    email = html.find('div', id="dzyj")
    if email:
        email = email.get_text().split('"')[1]
    else:
        email = '-'
    sleep(1)
    print '\t'.join([tutor_name, email, phone, url])
    return [tutor_name, email, phone, url]
#--------------------------------

def parseEachTutorIOZ(url):
    html = BeautifulSoup(br.open(url), 'html.parser')
    table = html.find('table', attrs={'width': '540px'})
    if not table:
        print >>sys.stderr, "No table found for "+ url
        return
    tr = table.find_all('tr')
    if not tr:
        print >>sys.stderr, "No tr found for "+ url
        return
    #for row in tr:
    #    tdL = row.find_all('td')
    tutor_name = tr[0].find_all('td')
    if tutor_name:
        tutor_name = tutor_name[1].get_text()
    else:
        tutor_name = '-'
    phone = tr[2].find_all('td')
    if phone:
        phone = phone[1].get_text()
    else:
        phone = '-'
    email = tr[3].find_all('td')
    if email:
        email = email[1].get_text()
    else:
        email = '-'
    sleep(1)
    print '\t'.join([tutor_name, email, phone, url])
    return [tutor_name, email, phone, url]
#--------------------------------

def parseEachTutorIB(url):
    html = BeautifulSoup(br.open(url), 'html.parser')
    table = html.find('table', attrs={'id': 'table24'})
    if not table:
        print >>sys.stderr, "No table found for "+ url
        return
    tr = table.find_all('tr')
    if not tr:
        print >>sys.stderr, "No tr found for "+ url
        return
    #for row in tr:
    #    tdL = row.find_all('td')
    tutor_name = tr[2].find_all('td')
    if tutor_name:
        tutor_name = tutor_name[0].get_text()
    else:
        tutor_name = '-'

    table = html.find('table', attrs={'id': 'table25'})
    if not table:
        print >>sys.stderr, "No table found for "+ url
        return
    tr = table.find_all('tr')
    if not tr:
        print >>sys.stderr, "No tr found for "+ url
        return
    phone = tr[1].find_all('td')
    if phone:
        phone = phone[0].get_text()
    else:
        phone = '-'
    email = tr[2].find_all('td')
    if email:
        email = email[0].find('a').get('href').replace('mailto:', '')
    else:
        email = '-'
    sleep(1)
    print '\t'.join([tutor_name, email, phone, url])
    return [tutor_name, email, phone, url]
#--------------------------------
def parseTsingHua(url):
    nite_A = br.open(url)
    html = BeautifulSoup(nite_A,  'html.parser')
    ul = html.find("ul", class_="list clearfloat")
    liL = ul.find_all("li")
    for li in liL:
        tutor_link = li.find('a').get('href')
        pL = li.find_all('p')
        tutor_name = pL[0].get_text().replace(', PhD', '')
        tutor_mail = pL[1].get_text().split('\n')[-1]
        print '\t'.join([tutor_name, tutor_mail, '-', tutor_link])
#--------------------------------------
def parseEachTutorTsingHuaLife(tutor_url):
    #http://life.tsinghua.edu.cn/faculty/faculty/329.html
    html = BeautifulSoup(br.open(tutor_url), 'html.parser')
    table = html.find('table', attrs={"width":"594"})
    tutor_name = table.find('td', class_="STYLE20").get_text()
    table = html.find('table', attrs={"width":"620"})
    table_name = table.find('td')
    try:
        table_name = table_name.find_all('div')[-1].get_text()
        phone, mail = table_name.split()
    except IndexError:
        table_name = table_name.find_all('p')
        phone = table_name[-2].get_text()
        mail = str(table_name[-1].get_text())
    except:
        print >>sys.stderr, tutor_url
        print >>sys.stderr, table
        #sys.exit(1)
        return
    print '\t'.join([tutor_name, mail, phone, tutor_url])


#-----------------------------------------

def parseTsingHuaLife(url):
    #http://life.tsinghua.edu.cn/faculty/faculty/229.html
    html = BeautifulSoup(br.open(url),  'html.parser')
    table = html.find('table', attrs={"width":"590"})
    aL = table.find_all('a')
    tutor_urlL = [i.get('href') for i in aL]

    tutor_urlL = set(tutor_urlL)
    tutor_infoL = []
    for tutor_url in tutor_urlL:
        tutor_infoL.append(parseEachTutorTsingHuaLife(tutor_url))
    return tutor_infoL
    
#-------------------------------------

def parseIGDB(url):
    nite_A = br.open(url)
    nite_A_html = BeautifulSoup(nite_A,  'html.parser')
    tutor_url  = nite_A_html.find_all('a', class_="t2_link")
    tutor_urlL = [i.get('href') for i in tutor_url]

    tutor_urlL = set(tutor_urlL)
    tutor_infoL = []
    for tutor_url in tutor_urlL:
        tutor_infoL.append(parseEachTutorIGDB(tutor_url))
    return tutor_infoL
#-----------------------------------------------

def parseCAS(url):
    nite_A = br.open(url)
    nite_A_html = BeautifulSoup(nite_A,  'html.parser')
    tutor_url  = nite_A_html.find_all('a')
    tutor_urlL = [i.get('href') for i in tutor_url]
    tutor_urlL = [i for i in tutor_urlL if i.find('sourcedb')!=-1]

    tutor_urlL = set(tutor_urlL)
    return tutor_urlL
#---------------------------------------

def parseIOZ(url):
    tutor_urlL = parseCAS(url)
    tutor_infoL = []
    for tutor_url in tutor_urlL:
        tutor_infoL.append(parseEachTutorIOZ(tutor_url))
    return tutor_infoL
#---------------------

def parseIB(url):
    tutor_urlL = parseCAS(url)
    tutor_infoL = []
    for tutor_url in tutor_urlL:
        tutor_infoL.append(parseEachTutorIB(tutor_url))
    return tutor_infoL
#---------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    fileL = [i.strip() for i in file.split(',')]
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    urlD = {'IGDB': 
               {"url": "http://www.genetics.ac.cn/rcjy/", 
                "func": parseIGDB}, 
            'IOZ':
                {"url": "http://www.ioz.cas.cn/rcjy/", 
                "func": parseIOZ
                }, 
            'IB':
                {"url": "http://www.ib.cas.cn/duiwu/zuzhang/", 
                "func": parseIB
                }, 
                'TsingHua': {"url":"http://www.sps.tsinghua.edu.cn/cn/team/team/", 
                    "func": parseTsingHua}, 
                'TsingHuaLife': {"url": "http://life.tsinghua.edu.cn/faculty/faculty/229.html", 
                    "func": parseTsingHuaLife
                    }, 
                'Daoshi_1':{'url': 'parse_cksp_eol_cn', 
                    "func": parse_cksp_eol_cn
                    }
            }

    for file in fileL:
        instituteD = urlD.get(file, '')
        if instituteD:
            tutor_infoL = instituteD['func'](instituteD['url'])
            #print "\n".join(['\t'.join(i) for i in tutor_infoL])

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


