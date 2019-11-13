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
    This is designed to detect the datatype of each column of data matrix and generate a SQL command for creating data tables.

    1. The first line will be treated as header line. All non-alphabetical words (except <_> or numbers) will be transfeered to '_'.
    2. If all items are strings in one column, the maximum length will be extracted.
    3. If all items are int or float, special numerical type will be generated.
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
import pandas as pd
import re

import numpy as np
import cStringIO
import pandas.io.sql as psql
from dateutil import parser

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
        metavar="FILEIN", help="Data matrix")
    parser.add_option("-n", "--data-table-name", dest="data_table_name",
        metavar="FILEIN", help="Name for data table")
    parser.add_option("-m", "--max-string_length", dest="max_str_len",
        type="int", default=21840, help="If length of string larger than given value (default 21840), this column will be treated as <text>.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

dbtypes={
    'mysql' : {'DATE':'DATE', 'DATETIME':'DATETIME',           'INT':'BIGINT',  'FLOAT':'FLOAT',  'VARCHAR':'VARCHAR'},
    'oracle': {'DATE':'DATE', 'DATETIME':'DATE',               'INT':'NUMBER',  'FLOAT':'NUMBER', 'VARCHAR':'VARCHAR2'},
    'sqlite': {'DATE':'TIMESTAMP', 'DATETIME':'TIMESTAMP',     'INT':'NUMBER',  'FLOAT':'NUMBER', 'VARCHAR':'VARCHAR2'},
    'postgresql': {'DATE':'TIMESTAMP', 'DATETIME':'TIMESTAMP', 'INT':'BIGINT',  'FLOAT':'REAL',   'VARCHAR':'TEXT'},
}

# from read_frame.  ?datetime objects returned?  convert to datetime64?
def read_db(sql, con):
    return []
    #return psql.frame_query(sql, con)


def table_exists(name=None, con=None, flavor='sqlite'):
    if flavor == 'sqlite':
        sql="SELECT name FROM sqlite_master WHERE type='table' AND name='MYTABLE';".replace('MYTABLE', name)
    elif flavor == 'mysql':
        sql="show tables like 'MYTABLE';".replace('MYTABLE', name)
    elif flavor == 'postgresql':
        sql= "SELECT * FROM pg_tables WHERE tablename='MYTABLE';".replace('MYTABLE', name)
    elif flavor == 'oracle':
        sql="select table_name from user_tables where table_name='MYTABLE'".replace('MYTABLE', name.upper())
    elif flavor == 'odbc':
        raise NotImplementedError
    else:
        raise NotImplementedError
    
    df = read_db(sql, con)
    print sql, df
    print 'table_exists?', len(df)
    exists = True if len(df)>0 else False
    return exists

def write_frame(frame, name=None, con=None, flavor='sqlite', if_exists='fail'):
    """
    Write records stored in a DataFrame to specified dbms. 
    
    if_exists:
        'fail'    - create table will be attempted and fail
        'replace' - if table with 'name' exists, it will be deleted        
        'append'  - assume table with correct schema exists and add data.  if no table or bad data, then fail.
            ??? if table doesn't exist, make it.
        if table already exists.  Add: if_exists=('replace','append','fail')
    """

    if if_exists=='replace' and table_exists(name, con, flavor):    
        cur = con.cursor()   
        cur.execute("drop table "+name)
        cur.close()    
    
    if if_exists in ('fail','replace') or ( if_exists=='append' and table_exists(name, con, flavor)==False ):
        #create table
        schema = get_schema(frame, name, flavor)
        if flavor=='oracle':
            schema = schema.replace(';','')
        cur = con.cursor()    
        if flavor=='mysql':
            cur.execute("SET sql_mode='ANSI_QUOTES';")
        print 'schema\n', schema
        cur.execute(schema)
        cur.close()
        print 'created table' 
        
    cur = con.cursor()
    #bulk insert
    if flavor=='sqlite' or flavor=='odbc':       
        wildcards = ','.join(['?'] * len(frame.columns))
        insert_sql = 'INSERT INTO %s VALUES (%s)' % (name, wildcards)
        #print 'insert_sql', insert_sql
        data = [tuple(x) for x in frame.values]
        #print 'data', data
        cur.executemany(insert_sql, data)
        
    elif flavor=='oracle':
        cols=[db_colname(k) for k in frame.dtypes.index]
        colnames = ','.join(cols)
        colpos = ', '.join([':'+str(i+1) for i,f in enumerate(cols)])
        insert_sql = 'INSERT INTO %s (%s) VALUES (%s)' % (name, colnames, colpos)
        #print 'insert_sql', insert_sql
        data = [ convertSequenceToDict(rec) for rec in frame.values] 
        #print data
        cur.executemany(insert_sql, data)
        
    elif flavor=='mysql':
        
        wildcards = ','.join(['%s'] * len(frame.columns))
        cols=[db_colname(k) for k in frame.dtypes.index]
        colnames = ','.join(cols)
        insert_sql = 'INSERT INTO %s (%s) VALUES (%s)' % (name, colnames, wildcards)
        print insert_sql
        #data = [tuple(x) for x in frame.values]
        data= [ tuple([ None if isnull(v) else v for v in rw]) for rw in frame.values ] 
        print data[0]
        cur.executemany(insert_sql, data)
        
    elif flavor=='postgresql':
        postgresql_copy_from(frame, name, con)    
    else:
        raise NotImplementedError        
    con.commit()
    cur.close()
    return

def nan2none(df):
    dnp = df.values
    for rw in dnp:
        rw2 = tuple([ None if v==np.Nan else v for v in rw])
        
    tpl_list= [ tuple([ None if v==np.Nan else v for v in rw]) for rw in dnp ] 
    return tpl_list
    
def db_colname(pandas_colname):
    '''convert pandas column name to a DBMS column name
        TODO: deal with name length restrictions, esp for Oracle
    '''
    colname =  pandas_colname.replace(' ','_').strip()                  
    return colname
    

def postgresql_copy_from(df, name, con ):
    # append data into existing postgresql table using COPY
    
    # 1. convert df to csv no header
    output = cStringIO.StringIO()
    
    # deal with datetime64 to_csv() bug
    have_datetime64 = False
    dtypes = df.dtypes
    for i, k in enumerate(dtypes.index):
        dt = dtypes[k]
        print 'dtype', dt, dt.itemsize
        if str(dt.type)=="<type 'numpy.datetime64'>":
            have_datetime64 = True

    if have_datetime64:
        d2=df.copy()    
        for i, k in enumerate(dtypes.index):
            dt = dtypes[k]
            if str(dt.type)=="<type 'numpy.datetime64'>":
                d2[k] = [ v.to_pydatetime() for v in d2[k] ]                
        #convert datetime64 to datetime
        #ddt= [v.to_pydatetime() for v in dd] #convert datetime64 to datetime
        d2.to_csv(output, sep='\t', header=False, index=False)
    else:
        df.to_csv(output, sep='\t', header=False, index=False)                        
    output.seek(0)
    contents = output.getvalue()
    print 'contents\n', contents
       
    # 2. copy from
    cur = con.cursor()
    cur.copy_from(output, name)    
    con.commit()
    cur.close()
    return


#source: http://www.gingerandjohn.com/archives/2004/02/26/cx_oracle-executemany-example/
def convertSequenceToDict(list):
    """for  cx_Oracle:
        For each element in the sequence, creates a dictionary item equal
        to the element and keyed by the position of the item in the list.
        >>> convertListToDict(("Matt", 1))
        {'1': 'Matt', '2': 1}
    """
    dict = {}
    argList = range(1,len(list)+1)
    for k,v in zip(argList, list):
        dict[str(k)] = v
    return dict

    
def get_schema(frame, name, flavor, varchar_max_len):
    types = dbtypes[flavor]  #deal with datatype differences
    column_types = []
    dtypes = frame.dtypes
    for i,k in enumerate(dtypes.index):
        dt = dtypes[k]
        #print 'dtype', dt, dt.itemsize
        if str(dt.type)=="<type 'numpy.datetime64'>":
            sqltype = types['DATETIME']
        elif issubclass(dt.type, np.datetime64):
            sqltype = types['DATETIME']
        elif issubclass(dt.type, (np.integer, np.bool_)):
            sqltype = types['INT']
        elif issubclass(dt.type, np.floating):
            sqltype = types['FLOAT']
        else:
            #sampl = frame[ frame.columns[i] ][0]
            sampl = frame.iloc[0,i]
            #print 'other', type(sampl)    
            if str(type(sampl))=="<type 'datetime.datetime'>":
                sqltype = types['DATETIME']
            elif str(type(sampl))=="<type 'datetime.date'>":
                sqltype = types['DATE']                   
            else:
                if flavor in ('mysql','oracle'):                
                    size = 2 + max( (len(str(a)) for a in frame[k]) )
                    print >>sys.stderr, k,'varchar sz', size
                    if size > varchar_max_len:
                        sqltype = 'TEXT'
                    else:
                        sqltype = types['VARCHAR'] + '(?)'.replace('?', str(size) )
                else:
                    sqltype = types['VARCHAR']
        colname =  db_colname(k)  #k.upper().replace(' ','_')                  
        column_types.append((colname, sqltype))
    columns = ',\n  '.join('%s %s' % x for x in column_types)
    template_create = """CREATE TABLE %(name)s (
    %(columns)s
    );"""    
    #print 'COLUMNS:\n', columns
    create = template_create % {'name' : name, 'columns' : columns}
    return create
    

###############################################################################

def test_sqlite(name, testdf):
    print '\nsqlite, using detect_types=sqlite3.PARSE_DECLTYPES for datetimes'
    import sqlite3
    with sqlite3.connect('test.db', detect_types=sqlite3.PARSE_DECLTYPES) as conn:
        #conn.row_factory = sqlite3.Row
        write_frame(testdf, name, con=conn, flavor='sqlite', if_exists='replace')
        df_sqlite = read_db('select * from '+name, con=conn)    
        print 'loaded dataframe from sqlite', len(df_sqlite)   
    print 'done with sqlite'


def test_oracle(name, testdf):
    print '\nOracle'
    import cx_Oracle
    with cx_Oracle.connect('YOURCONNECTION') as ora_conn:
        testdf['d64'] = np.datetime64( testdf['hire_date'] )
        write_frame(testdf, name, con=ora_conn, flavor='oracle', if_exists='replace')    
        df_ora2 = read_db('select * from '+name, con=ora_conn)    

    print 'done with oracle'
    return df_ora2
   
    
def test_postgresql(name, testdf):
    #from pg8000 import DBAPI as pg
    import psycopg2 as pg
    print '\nPostgresQL, Greenplum'    
    pgcn = pg.connect(YOURCONNECTION)
    print 'df frame_query'
    try:
        write_frame(testdf, name, con=pgcn, flavor='postgresql', if_exists='replace')   
        print 'pg copy_from'    
        postgresql_copy_from(testdf, name, con=pgcn)    
        df_gp = read_db('select * from '+name, con=pgcn)    
        print 'loaded dataframe from greenplum', len(df_gp)
    finally:
        pgcn.commit()
        pgcn.close()
    print 'done with greenplum'

 
def test_mysql(name, testdf):
    import MySQLdb
    print '\nmysql'
    cn= MySQLdb.connect(YOURCONNECTION)    
    try:
        write_frame(testdf, name='test_df', con=cn, flavor='mysql', if_exists='replace')
        df_mysql = read_db('select * from '+name, con=cn)    
        print 'loaded dataframe from mysql', len(df_mysql)
    finally:
        cn.close()
    print 'mysql done'
#------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    max_str_len = options.max_str_len
    data_table_name = options.data_table_name
    database_type = 'mysql'
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    #data = pd.read_table(file, sep="\t" , header=0, index_col=0)
    data = pd.read_table(file, sep="\t" , header=0, index_col=0)
    colnames_ori = data.columns
    newCol_name = [re.sub('[^A-Za-z_]', '_', i) for i in colnames_ori]
    data.columns = newCol_name
    print get_schema(data, data_table_name, database_type, max_str_len)
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


