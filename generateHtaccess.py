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
    This is designed to generate htpasswd and .htaccess file.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

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
    parser.add_option("-u", "--user", dest="user",
        metavar="USER-NAME", help="Username")
    parser.add_option("-p", "--passwd", dest="passwd",
        metavar="PASSWORD", help="PASSWORD")
    parser.add_option("-d", "--secret-file-dir-local", dest="secrect_file_dir_local",
        help="Directory for file to save username and password info.")
    parser.add_option("-s", "--secret-file-name", dest="secrect_file_name",
        help="Name for file to save username and password info.")
    parser.add_option("-D", "--secret-file-dir-server", dest="secrect_file_dir_server",
        help="Directory for file to save username and password info at remote server.")
    parser.add_option("-i", "--ip", dest="ip",
        help="< >separated list of IPs have full authority")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    (options, args) = parser.parse_args(argv[1:])
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    user = options.user
    passwd = options.passwd
    secrect_file_name = options.secrect_file_name
    secrect_file_dir_local = options.secrect_file_dir_local.rstrip('/')+'/'
    secrect_file_dir_server = options.secrect_file_dir_server.rstrip('/')+'/'
    ip = options.ip
    cmd = ["mkdir -p", secrect_file_dir_local]
    os.system(' '.join(cmd))
    cmd = ['htpasswd -cmb', secrect_file_dir_local+secrect_file_name, user, passwd]
    os.system(' '.join(cmd))

    print '''<If "%{{HTTP_User_Agent}} !~ /MicroMessenger/">
AuthType Basic
AuthName 'Please contact ct@ehbio.com for accssing.'
AuthUserFile {}{}
Require valid-user
</If>   
'''.format(secrect_file_dir_server, secrect_file_name)
#----------------------------------------

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


