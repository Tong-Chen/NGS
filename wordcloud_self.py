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
    This is designed to generate word cloud for WeChat.
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
#reload(sys)
#sys.setdefaultencoding('utf8')

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
        metavar="FILEIN", help="A text file containing blog text")
    parser.add_option("-s", "--stopwords", dest="stopwords",
        default="/home/ct/blog/tong-chen.github.com/_posts/stopwords",
        help="A file containing words to be ignored. One each row. Default system default.")
    parser.add_option("-p", "--pic", dest="pic",
        default="~/blog/tong-chen.github.com/_posts/logo_mode.png", 
        help="A picture to be moded on. Default system default.")
    parser.add_option("-w", "--words-extra", dest="words",
            help="SPecial words to be added to input file to weight it. In format like <airflow:30;docker:40> meaning add 30 airflow words and 40 docker words.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    stopwords = options.stopwords
    pic = options.pic
    words = options.words
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    output_pic = file + '.png'
    if words:
        wordsL = [i.split(':') for i in words.split(';')]
        output_tmp = file + '.tmp_wordcloud_747474'
        fh_out = open(output_tmp, 'w')
        #--------------------------------
        for line in open(file):
            print >>fh_out, line,
        #-------------END reading file----------
        for word, cnt in wordsL:
            for i in range(int(cnt)):
                print >>fh_out, word
        fh_out.close()
        #----close file handle for files-----
    else:
        output_tmp = file
    cmd = ["wordcloud_cli.py", "--text", output_tmp, "--imagefile", output_pic, "--stopwords",  stopwords, "--mask", pic, "--colormask", pic]
    os.system(' '.join(cmd))
    if output_tmp != file:
        os.system("/bin/rm -f "+output_tmp)
    #-----------end close fh-----------
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


