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
    This will generate a dag file for uploading data to GEO.
    
    Suppose we have data for uploading in directory `/home/name/project1/upload`, 
    the generated dag file would be saved at the directory `/home/name/project1/`.
'''

import sys
import os
from json import dumps as json_dumps
from datetime import datetime, timedelta 
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
    usages = "%prog -i ftp"
    parser = OP(usage=usages)
    parser.add_option("-f", "--ftp-link", dest="ftp",
        metavar="FTP", default="ftp-private.ncbi.nlm.nih.gov", 
        help="The url of FTP server. Default <ftp-private.ncbi.nlm.nih.gov>.")
    parser.add_option("-u", "--username-ftp", dest="userFTP",
        metavar="USERNAME", default='geo', 
        help="Username of FTP server. Default <geo>.")
    parser.add_option("-p", "--passwd", dest="passwd",
        metavar="PASSWORD", default='33%9uyj_fCh?M16H', 
        help="Password for FTP server. Default <33%9uyj_fCh?M16H>.")
    parser.add_option("-U", "--user-name-NCBI", dest="userNCBI",
        help="Username for NCBI login.")
    parser.add_option("-d", "--absolute-path-for-directory-containing-all-data", 
        dest="dir",  
        help="The absolute path for directory containing \
all and only data for uploading. \n\
Please make sure you have another copy of all data in this folder, \
since all these data will be deleted after uploading.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    (options, args) = parser.parse_args(argv[1:])
    assert options.ftp != None, "A ftp link needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    ftp = options.ftp
    userFTP = options.userFTP
    passwd  = options.passwd
    userNCBI = options.userNCBI
    verbose = options.verbose
    dir = options.dir
    parentdir = dir.rstrip('/').rsplit('/', 1)[0]
    output_dag = open(parentdir+'/'+userNCBI + '.dag.py', 'w')
    lftp_script1 = parentdir+'/'+userNCBI + '_lftp.script1'
    lftp_script2 = parentdir+'/'+userNCBI + '_lftp.script2'
    output_lftp1 = open(lftp_script1, 'w')
    output_lftp2 = open(lftp_script2, 'w')
    
    global debug
    #-----------------------------------
    print >>output_lftp1, '''
open -u %s,%s %s
mkdir -p fasp/GEO_metadata_%s/
    ''' % (userFTP, passwd, ftp, userNCBI)

    print >>output_lftp2, '''
open -u %s,%s %s
cd fasp/GEO_metadata_%s/
cache size 33554432
set cmd:parallel 20
mput -E %s/*
    ''' % (userFTP, passwd, ftp, userNCBI, dir)

    output_lftp1.close()
    output_lftp2.close()

    print >>output_dag, '''
from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from airflow.operators.email_operator import EmailOperator


from datetime import datetime, timedelta


one_min_ago = datetime.combine(datetime.today() - timedelta(minutes=20),
                                  datetime.min.time())

default_args = {
    'owner': 'ct',         
    'depends_on_past': True, 
    'start_date': one_min_ago,
    'email': ['chentong_biology@163.com'],
    'email_on_failure': True, 
    'email_on_retry': True, 
    'retries': 500, 
    'retry_delay': timedelta(seconds=30)
}
'''

    print >>output_dag, '''

dag_id = "%s_GEO_upload"

mkdir = "lftp -f %s"
lftp = "lftp -f %s"


dag = DAG(dag_id, default_args=default_args,
    schedule_interval="@once")

geo_lftp_mkdir = BashOperator(
    task_id='geo_lftp_mkdir', 
    bash_command=mkdir, 
    dag=dag, retries=3)

#In case the directory exists, but still supplying 3 minutes for `geo_lftp_mkdir` to run before strating `geo_lftp_upload`
geo_lftp_sleep = BashOperator(
    task_id='geo_lftp_sleep', 
    bash_command="sleep 3m",
    dag=dag, retries=1, trigger_rule="dummy")

geo_lftp_sleep.set_upstream(geo_lftp_mkdir)

geo_lftp_upload = BashOperator(
    task_id='geo_lftp_upload', 
    bash_command=lftp, 
    dag=dag)

geo_lftp_upload.set_upstream(geo_lftp_sleep)

''' % (userNCBI, lftp_script1, lftp_script2)

    twoyearafter = datetime.now()+timedelta(days=750)
    
    print >>output_dag, """
html_content = '''
Dear Sir/Madam,  
  
  Thanks for you kindly host such great public data resource.
    
  I have successfully transferred my data to NCBI-GEO ftp sever. 
      
  Here is the information you may be needed for further processing
        
    1. GEO account username: <strong>%s</strong>
    2. Names of the directory and files deposited: <strong>fasp/GEO_metadata_%s/</strong>
    3. Public release date: <strong>%s</strong>
      
  If there is any format or content problem, please do not hesitate to contact me.
            
  Best,  
              
  %s

''' 
""" % (userNCBI, userNCBI, twoyearafter.strftime('%Y-%m-%d'), userNCBI)
    print >>output_dag, '''

Success_mail = EmailOperator(
    task_id="Success_mail", 
    to="chentong_biology@163.com",  
    subject="ftp upload",  
    html_content=html_content,  
    dag=dag)

Success_mail.set_upstream(geo_lftp_upload)

'''

    output_dag.close()

if __name__ == '__main__':
    main()


