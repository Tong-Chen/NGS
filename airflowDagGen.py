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
    This is designed to generate airflow dag for batch assignments.

config.json

{
    "Attention": [
    "1. No comma sign for the last item in each group.",
    "2. No single quote or parenthesis is allowed.", 
    "3. No chinese comma out of double quote.",
    "4. Use list to order subjects."
    ],

    "work_dir": "/working-directory", 
    "program": "virtualScreening.py", 
    "program_parameter": [
        ["-p", ["prot1.pdb", "prot2.pdb", "prot3.pdb"]], 
        ["-l", ["chem1.pdb", "chem2.pdb"]], 
        ["-o", "result"]
    ], 

    "airflow_parameter":{
        "owner": "'ct'", 
        "dag_id" : "'vs'", 
        "start_date": "one_min_ago", 
        "depends_on_past": "False", 
        "email": ["chentong_biology@163.com"], 
        "email_on_failure": "True", 
        "email_on_retry": "True", 
        "retries": 500, 
        "retry_delay": "timedelta(hours=30)", 
        "schedule_interval": "'@once'",
        "comment": {
            "comment": "all airflow supported parameters.", 
            "start_date": "can be <one_min_ago> or <'2016-01-01'>"
            "special_attention": "sting value shuould be quoted in ***single quote (')***"
        }
    }, 

    "parallel_parameter": {
        "maxrun": 5, 
        "comment": {
            "maxrun": "Set the number of assignments running parallely"
        }
    }, 

    "comment": {
        "comment": "Explanation words (optional)", 
        "program": "The program wants to run", 
        "program_parameter": "Parameters given to the program. For parameters accepting one or more files, both a list as exampled above for <-l> or a string like <\"'chem1.pdb','chemb.pdb'\"> are acceptable.", 
        "work_dir": "The absolute path for working directory"
    }
}

generated DAG


from airflow import DAG
from airflow.operators import BashOperator, EmailOperator

from datetime import datetime, timedelta


one_min_ago = datetime.combine(datetime.today() + timedelta(minutes=20),
                                  datetime.min.time())

default_args = {
    'email_on_retry': True,
    'email': [u'chentong_biology@163.com'],
    'email_on_failure': True,
    'retry_delay': timedelta(hours=30),
    'owner': 'ct',
    'depends_on_past': False,
    'start_date': one_min_ago,
    'retries': 500
}


dag = DAG('vs', default_args=default_args, schedule_interval='@once')


chem1_pdb_prot1_pdb = BashOperator(
    task_id='chem1_pdb_prot1_pdb', 
    bash_command="(cd /working-directory; virtualScreening.py -l chem1.pdb -o result -p prot1.pdb) ", 
    dag=dag)

chem1_pdb_prot1_pdb_success_mail = EmailOperator(
    task_id="chem1_pdb_prot1_pdb_success_mail", 
    to=[u'chentong_biology@163.com'],  
    subject="chem1_pdb_prot1_pdb success",  
    html_content="chem1_pdb_prot1_pdb success",  
    dag=dag)
                
chem1_pdb_prot1_pdb_success_mail.set_upstream(chem1_pdb_prot1_pdb)
#chem1_pdb_prot1_pdb.set_upstream( )


chem1_pdb_prot2_pdb = BashOperator(
    task_id='chem1_pdb_prot2_pdb', 
    bash_command="(cd /working-directory; virtualScreening.py -l chem1.pdb -o result -p prot2.pdb) ", 
    dag=dag)

chem1_pdb_prot2_pdb_success_mail = EmailOperator(
    task_id="chem1_pdb_prot2_pdb_success_mail", 
    to=[u'chentong_biology@163.com'],  
    subject="chem1_pdb_prot2_pdb success",  
    html_content="chem1_pdb_prot2_pdb success",  
    dag=dag)
                
chem1_pdb_prot2_pdb_success_mail.set_upstream(chem1_pdb_prot2_pdb)
#chem1_pdb_prot2_pdb.set_upstream( )


chem1_pdb_prot3_pdb = BashOperator(
    task_id='chem1_pdb_prot3_pdb', 
    bash_command="(cd /working-directory; virtualScreening.py -l chem1.pdb -o result -p prot3.pdb) ", 
    dag=dag)

chem1_pdb_prot3_pdb_success_mail = EmailOperator(
    task_id="chem1_pdb_prot3_pdb_success_mail", 
    to=[u'chentong_biology@163.com'],  
    subject="chem1_pdb_prot3_pdb success",  
    html_content="chem1_pdb_prot3_pdb success",  
    dag=dag)
                
chem1_pdb_prot3_pdb_success_mail.set_upstream(chem1_pdb_prot3_pdb)
#chem1_pdb_prot3_pdb.set_upstream( )


chem2_pdb_prot1_pdb = BashOperator(
    task_id='chem2_pdb_prot1_pdb', 
    bash_command="(cd /working-directory; virtualScreening.py -l chem2.pdb -o result -p prot1.pdb) ", 
    dag=dag)

chem2_pdb_prot1_pdb_success_mail = EmailOperator(
    task_id="chem2_pdb_prot1_pdb_success_mail", 
    to=[u'chentong_biology@163.com'],  
    subject="chem2_pdb_prot1_pdb success",  
    html_content="chem2_pdb_prot1_pdb success",  
    dag=dag)
                
chem2_pdb_prot1_pdb_success_mail.set_upstream(chem2_pdb_prot1_pdb)
#chem2_pdb_prot1_pdb.set_upstream( )


chem2_pdb_prot2_pdb = BashOperator(
    task_id='chem2_pdb_prot2_pdb', 
    bash_command="(cd /working-directory; virtualScreening.py -l chem2.pdb -o result -p prot2.pdb) ", 
    dag=dag)

chem2_pdb_prot2_pdb_success_mail = EmailOperator(
    task_id="chem2_pdb_prot2_pdb_success_mail", 
    to=[u'chentong_biology@163.com'],  
    subject="chem2_pdb_prot2_pdb success",  
    html_content="chem2_pdb_prot2_pdb success",  
    dag=dag)
                
chem2_pdb_prot2_pdb_success_mail.set_upstream(chem2_pdb_prot2_pdb)
#chem2_pdb_prot2_pdb.set_upstream( )


chem2_pdb_prot3_pdb = BashOperator(
    task_id='chem2_pdb_prot3_pdb', 
    bash_command="(cd /working-directory; virtualScreening.py -l chem2.pdb -o result -p prot3.pdb) ", 
    dag=dag)

chem2_pdb_prot3_pdb_success_mail = EmailOperator(
    task_id="chem2_pdb_prot3_pdb_success_mail", 
    to=[u'chentong_biology@163.com'],  
    subject="chem2_pdb_prot3_pdb success",  
    html_content="chem2_pdb_prot3_pdb success",  
    dag=dag)
                
chem2_pdb_prot3_pdb_success_mail.set_upstream(chem2_pdb_prot3_pdb)
chem2_pdb_prot3_pdb.set_upstream(chem1_pdb_prot1_pdb)



'''

import sys
import os
from json import dumps as json_dumps
from json import loads as json_loads
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
from string import maketrans
import re

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
        metavar="FILEIN", help="A config file in JSON format as described above.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def iterateParameters(optL, cmdD):
    '''
    optL = [
            ["-p prot1", "-p prot2", "-p prot3"], 
            ["-l chem1", "-l chem2"], 
            ["-o result"]
           ]

    cmdD = {
        "-p prot1": {"-l chem1": "result", "-l chem2": "result"}
    }

    '''
    if debug:
        print >>sys.stderr, "New"
        print >>sys.stderr, optL
        print >>sys.stderr, cmdD
    for opt in optL[0]:
        newOptL = optL[1:]
        if len(newOptL)>1:
            tmpD = {}
            cmdD[opt] = tmpD
            if debug:
                print >>sys.stderr, "Next"
                print >>sys.stderr, newOptL
                print >>sys.stderr, cmdD
            iterateParameters(newOptL, tmpD)
        else:
            cmdD[opt] = {}
            for last_opt in newOptL[0]:
                cmdD[opt][last_opt] = ''
#------------------------------------

def iterateParameters2(optL, recordL, cmdL):
    '''
    optL = [
            ["-p prot1", "-p prot2", "-p prot3"], 
            ["-l chem1", "-l chem2"], 
            ["-o result"]
           ]
    
    recordL: tmp record variable

    cmdL = [
        ["-p prot1", "-l chem1",  "result"], 
        ["-p prot1", "-l chem2",  "result"]
    ]

    '''
    if debug:
        print >>sys.stderr, "New"
        print >>sys.stderr, optL
        print >>sys.stderr, recordL
    for opt in optL[0]:
        newOptL = optL[1:]
        if len(newOptL) > 1:
            tmpL = recordL[:]
            tmpL.append(opt)
            if debug:
                print >>sys.stderr, "Next"
                print >>sys.stderr, newOptL
                print >>sys.stderr, tmpL
            iterateParameters2(newOptL, tmpL, cmdL)
        else:
            for last_opt in newOptL[0]:
                tmpL = recordL[:]
                tmpL.append(opt)
                tmpL.append(last_opt)
                cmdL.append(tmpL)
#------------------------------------


def generateCmd(program, program_parameterD, work_dir, maxrun, email):
    count = 1
    chdir = "(cd {}; ".format(work_dir)
    assignmentL = [' '] * maxrun

    optL = [[program]]
    nameL = []
    for opt, optV in program_parameterD.items():
        cmdL = []
        subNameL = []
        if type(optV) == list:
            for eachV in optV:
                cmd = ' '.join([opt, str(eachV)])
                cmdL.append(cmd)
                subNameL.append(str(eachV))
        else:
            cmdL = [' '.join([opt, str(optV)])]
        optL.append(cmdL)
        if subNameL:
            nameL.append(subNameL)
    #=------------------------------------------
    if debug:
        print >>sys.stderr, optL
        print >>sys.stderr, nameL
    recordL = []
    cmdLL = []
    iterateParameters2(optL, recordL, cmdLL)
    nameLL = []
    recordL = []
    iterateParameters2(nameL, recordL, nameLL)

    if debug:
        print >>sys.stderr, cmdLL
        print >>sys.stderr, nameLL

    num_start = re.compile('^[0-9]')
    illegal_char = re.compile('[^a-zA-Z0-9]')
    count = 1
    for cmdL in cmdLL:
        cmd = chdir+' '.join(cmdL)+')'
        task_name = "_".join(nameLL[count-1])
        task_name = illegal_char.sub('_', task_name)
        if num_start.match(task_name):
            task_name = '_'+task_name
        comment = '' if count > maxrun else '#'
        index = count%maxrun-1
        upstream_task = assignmentL[index]
        assignmentL[index] = task_name
        tmpD = {"cmd":cmd, "task_name":task_name, 
                "upstream_task":upstream_task, 
                "email": email, 
                "comment": comment}
        count += 1
        print '''
{d[task_name]} = BashOperator(
    task_id='{d[task_name]}', 
    bash_command="{d[cmd]} ", 
    dag=dag)

{d[task_name]}_success_mail = EmailOperator(
    task_id="{d[task_name]}_success_mail", 
    to={d[email]},  
    subject="{d[task_name]} success",  
    html_content="{d[task_name]} success",  
    dag=dag)
                
{d[task_name]}_success_mail.set_upstream({d[task_name]})
{d[comment]}{d[task_name]}.set_upstream({d[upstream_task]})
'''.format(d=tmpD)

#---------END of generateCmd----------------------

def generateCmd2(program, program_parameterL, work_dir, maxrun, email):
    '''
    program_parameterL = [
        ['-p', ['1.pdb', '2.pdb', '3.pdb']],
        ['-l', ['c1.pdb', 'c2.pdb']], 
        ['-o', 'result']
    ]
    '''
    count = 1
    chdir = "(cd {}; ".format(work_dir)
    assignmentL = [' '] * maxrun

    optL = [[program]]
    nameL = []
    for opt, optV in program_parameterL:
        cmdL = []
        subNameL = []
        if type(optV) == list:
            for eachV in optV:
                cmd = ' '.join([opt, str(eachV)])
                cmdL.append(cmd)
                subNameL.append(str(eachV))
        else:
            cmdL = [' '.join([opt, str(optV)])]
        optL.append(cmdL)
        if subNameL:
            nameL.append(subNameL)
    #=------------------------------------------
    if debug:
        print >>sys.stderr, optL
        print >>sys.stderr, nameL
    recordL = []
    cmdLL = []
    iterateParameters2(optL, recordL, cmdLL)
    nameLL = []
    recordL = []
    iterateParameters2(nameL, recordL, nameLL)

    if debug:
        print >>sys.stderr, cmdLL
        print >>sys.stderr, nameLL

    num_start = re.compile('^[0-9]')
    illegal_char = re.compile('[^a-zA-Z0-9]')
    count = 1
    for cmdL in cmdLL:
        cmd = chdir+' '.join(cmdL)+')'
        task_name = "_".join(nameLL[count-1])
        task_name = illegal_char.sub('_', task_name)
        if num_start.match(task_name):
            task_name = '_'+task_name
        comment = '' if count > maxrun else '#'
        index = count%maxrun-1
        upstream_task = assignmentL[index]
        assignmentL[index] = task_name
        tmpD = {"cmd":cmd, "task_name":task_name, 
                "upstream_task":upstream_task, 
                "email": email, 
                "comment": comment}
        count += 1
        print '''
{d[task_name]} = BashOperator(
    task_id='{d[task_name]}', 
    bash_command="{d[cmd]} ", 
    dag=dag)

{d[task_name]}_success_mail = EmailOperator(
    task_id="{d[task_name]}_success_mail", 
    to={d[email]},  
    subject="{d[task_name]} success",  
    html_content="{d[task_name]} success",  
    dag=dag)
                
{d[task_name]}_success_mail.set_upstream({d[task_name]})
{d[comment]}{d[task_name]}.set_upstream({d[upstream_task]})
'''.format(d=tmpD)

#---------END of generateCmd----------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    config = options.filein
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    fh = open(config)
    context = ''.join(fh.readlines())
    fh.close()
    if debug:
        print >>sys.stderr, context
    configD = json_loads(context)
    #--------------------------------
    print '''
from airflow import DAG
from airflow.operators import BashOperator, EmailOperator

from datetime import datetime, timedelta


one_min_ago = datetime.combine(datetime.today() - timedelta(minutes=20),
                                  datetime.min.time())


default_args = {'''

    airflowArgD = configD['airflow_parameter']
    arg_keyL = airflowArgD.keys()
    arg_keyL.remove('dag_id')
    arg_keyL.remove('schedule_interval')
    if 'comment' in arg_keyL: 
        arg_keyL.remove('comment')
    len_argKey = len(arg_keyL)
    for i in range(1, len_argKey):
        key   = arg_keyL[i]
        value = airflowArgD[key]
        print "    '{}': {},".format(key, value)

    print "    '{}': {}\n{}".format(
        arg_keyL[0], airflowArgD[arg_keyL[0]], '}')
    print 

    print '''
dag = DAG({}, default_args=default_args, schedule_interval={})
'''.format(airflowArgD['dag_id'], airflowArgD['schedule_interval'])

    maxrun = configD["parallel_parameter"]["maxrun"]
    program = configD['program']
    program_parameter = configD['program_parameter']
    work_dir = configD['work_dir']
    generateCmd2(program, program_parameter, work_dir, maxrun, email=airflowArgD["email"])



if __name__ == '__main__':
    main()
