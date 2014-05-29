#!/usr/bin/python
# -*- coding: utf-8 -*-  

desc='''
Send emails.
'''

import sys
import os
from optparse import OptionParser as OP

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-r", "--receiver", dest='receiver',
            metavar="RECEIVER", default="chentong_biology@163.com", 
            help="The receiver for your mail. \
Default <chentong_biology@163.com>.")
    parser.add_option("-s", "--subject", dest='subject',
            metavar="SUBJECT", help="The subject for your mail")
    parser.add_option("-c", "--context", dest='context',
            metavar="CONTEXT", help="The main context for your mail. \
Symbol <-> should be given if you want SYS.STDIN as the mail content.")
    parser.add_option("-a", "--attachement", dest='attachment',
            metavar="ATTACHEMENT", help="The attachment for your mail. \
Multiple attachments should be separated by ','.")
    (options, args) = parser.parse_args(argv[1:])
    return (options, args)


import smtplib  
import mimetypes  
from email.MIMEMultipart import MIMEMultipart  
from email.MIMEBase import MIMEBase  
from email.MIMEText import MIMEText  
from email.MIMEAudio import MIMEAudio  
from email.MIMEImage import MIMEImage  
from email.Encoders import encode_base64  
  
def sendMail(subject, text, recipient, attachmentFilePaths):  
    gmailUser = 'mercury_gold@163.com'  
    gmailPassword = 'mercury.gold'  
  
    msg = MIMEMultipart()  
    msg['From'] = gmailUser  
    msg['To'] = recipient  
    msg['Subject'] = subject  
    msg.attach(MIMEText(text))  
  
    for attachmentFilePath in attachmentFilePaths:  
        if attachmentFilePath:
            msg.attach(getAttachment(attachmentFilePath))  
  
    mailServer = smtplib.SMTP('smtp.163.com')  
    mailServer.ehlo()  
    mailServer.starttls()  
    mailServer.ehlo()  
    mailServer.login(gmailUser, gmailPassword)  
    mailServer.sendmail(gmailUser, recipient, msg.as_string())  
    mailServer.close()  
  
    print('Sent email to %s' % recipient)  
  
def getAttachment(attachmentFilePath):  
    contentType, encoding = mimetypes.guess_type(attachmentFilePath)  
  
    if contentType is None or encoding is not None:  
        contentType = 'application/octet-stream'  
  
    mainType, subType = contentType.split('/', 1)  
    file = open(attachmentFilePath, 'rb')  
  
    if mainType == 'text':  
        attachment = MIMEText(file.read())  
    elif mainType == 'message':  
        attachment = email.message_from_file(file)  
    elif mainType == 'image':  
        attachment = MIMEImage(file.read(),_subType=subType)  
    elif mainType == 'audio':  
        attachment = MIMEAudio(file.read(),_subType=subType)  
    else:  
        attachment = MIMEBase(mainType, subType)  
    attachment.set_payload(file.read())  
    encode_base64(attachment)  
  
    file.close()  
  
    attachment.add_header('Content-Disposition', 'attachment',   filename=os.path.basename(attachmentFilePath))  
    return attachment  
  
  
# start to test  
def main():
    options, args = cmdparameter(sys.argv)
    subject = options.subject
    context = options.context
    if context != '-':
        content = context
    else:
        #fh = sys.stdin
        content = sys.stdin.read()
    recipient = options.receiver
    attach = [i.strip() for i in options.attachment.split(',')]

    sendMail(subject, content, recipient, attach) 

main()
