
from airflow import DAG
from airflow.operators import BashOperator, EmailOperator

from datetime import datetime, timedelta


one_min_ago = datetime.combine(datetime.today() - timedelta(minutes=20),
                                  datetime.min.time())

default_args = {
    'email_on_retry': True,
    'email': [u'chentong_biology@163.com'],
    'email_on_failure': True,
    'retry_delay': timedelta(seconds=30),
    'owner': 'ct',
    'depends_on_past': True,
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

