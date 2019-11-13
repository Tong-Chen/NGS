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
    Airflow based bat-running of I_TASSER.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
from string import maketrans
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
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="Fasta file with multiple protein sequences. Sequence name should not contain symbols other than alphabets, numbers and underlines.")
    parser.add_option("-m", "--maximum-number-itasser", dest="max_num",        
            type="int", help="Maximum allowed number of proteins for running I_TASSER paraellely.")
    parser.add_option("-n", "--project-name", dest="project_name",
            help="A string containing only alphabets, numbers and underlines to specify the name of the project. Normally one should add his/her username to avoid name conflicts.")
    parser.add_option("-w", "--working-dir", dest="work_dir",
        help="The absolute path for FASTA file")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def generateCmd(key, tmpL, count, max_num, index, assignmentL, work_dir):
    os.system("mkdir -p %s" % key)
    fh_out = open("%s/seq.fasta" % key, 'w')
    print >>fh_out, ">%s\n%s" % (key, ''.join(tmpL).rstrip('*'))
    fh_out.close()
    command = "(cd %s; I_TASSER.sh -n %s -d %s/)" % (work_dir, key, key)
    comment = '' if count > max_num else '#'

    print '''
%s = BashOperator(
    task_id='%s', 
    bash_command="%s", 
    dag=dag)

%s_success_mail = EmailOperator(
    task_id="%s_success_mail", 
    to="chentong_biology@163.com",  
    subject="%s I_TASSER success",  
    html_content="%s I_TASSER success",  
    dag=dag)
                
%s_success_mail.set_upstream(%s)
%s%s.set_upstream(%s)
''' % (key, key, command, 
       key, key, key, key, 
       key, key, 
       comment, key, assignmentL[index])
    if debug:
        print >>sys.stderr, "\t%s-->%s" % (key, assignmentL[index])
    assignmentL[index] = key
#---------END of generateCmd----------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    max_num = options.max_num
    name    = options.project_name
    work_dir = options.work_dir
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    print '''
from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from airflow.operators.email_operator import EmailOperator

from datetime import datetime, timedelta


one_min_ago = datetime.combine(datetime.today() - timedelta(minutes=20),
                                  datetime.min.time())

default_args = {
    'owner': 'ct',         
    'depends_on_past': False, 
    'start_date': one_min_ago,
    'email': ['chentong_biology@163.com'],
    'email_on_failure': True, 
    'email_on_retry': True, 
    'retries': 500, 
    'retry_delay': timedelta(hours=300)
}

dag = DAG('%s', default_args=default_args,
    schedule_interval="@once")

    ''' % (name)

    count = 1
    key = ''
    tmpL = ''
    assignmentL = ['_'] * max_num
    for line in fh:
        if line[0] == '>':
            if key:
                generateCmd(key, tmpL, count, max_num, index, assignmentL, work_dir)
                count += 1
                if debug:
                    print >>sys.stderr, assignmentL
            #---------------------------------------------------
            key = line[1:].strip().translate(maketrans(' .-/({[]})', '_'*10))
            key = '_'+key
            tmpL = []
            index = count % max_num - 1
            #print >>sys.stderr, index
        else:
            tmpL.append(line.strip())
    if key:
        generateCmd(key, tmpL, count, max_num, index, assignmentL, work_dir)
    #-------------END reading file----------
    #----close file handle for files-----
    if file != '-':
        fh.close()
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


