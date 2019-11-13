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
    This is designed to get the latest 10 articles of WeChat service platform by supplying WeChatID or key words.


logging.conf:

[loggers]
keys=root

[handlers]
keys=rotateFileHandler

[formatters]
keys=simpleFormatter

[logger_root]
level=WARNING
handlers=rotateFileHandler
qualname=simpleExample
propagate=0

[handler_rotateFileHandler]
class=handlers.RotatingFileHandler
level=WARNING
formatter=simpleFormatter
args=('log.txt', 'a+', 200000, 9)

[formatter_simpleFormatter]
format=%(asctime)s - [%(filename)s:%(lineno)d] - %(levelname)s - %(message)s
datefmt=

'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
reload(sys)
sys.setdefaultencoding('utf8')

import logging
import logging.config
from wechatsogou import *

debug = 0


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
        metavar="FILEIN", help="<;> separated multiple WeChatIDs. Default <jiweitianxia;biotrainee;BioGossip;bio_sxy;biobabble;xiaoddrz;hyyfreescience;genemeet;sw315ddayq;gene_test2016;>")
    parser.add_option("-k", "--keywords", dest="keywords",
        help="Multiple words separated by ';' for searching")
    parser.add_option("-l", "--loggingconf", dest="log_conf",
        default="/MPATHB/self/NGS/logging.conf", 
        help="Log file format. Default </MPATHB/self/NGS/logging.conf>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None || options.keywords != None, "Parameters needed for -i or -k"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    wechatid = options.filein

    if wechatid:
        wechatidL = wechatid.split(';')
    
    keywords = options.keywords

    if keywords:
        keywordsL = keywords.split(';')
    log_conf = options.log_conf
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    logging.config.fileConfig(log_conf)
    logger = logging.getLogger()
    wechats = WechatSogouApi()

    if wechatid:
        for wechat_id in wechatidL:
            #wechat_info = wechats.get_gzh_info(wechat_id)
            wechat_info = wechats.get_gzh_message(wechatid=wechat_id)
            

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


