import os
import sys
import time
import re

def copy(dir, *fileL):
    for file in fileL:
        if not os.path.exists(file):
            fh = open(file, 'w')
            print >>fh, '\t'.join(['No result' for i in range(10)])
            print >>fh, '\t'.join(['No result' for i in range(10)])
            fh.close()
            print >>sys.stderr, "** %s does not exist. We have written 10 columns with content <No result> to this file." % file
        #if os.path.exists(file):
        os.system("cp -u %s %s" % (file, dir))
#--------------------------------------------------------------

def copyzip(dir, *fileL):
    for file in fileL:
        if os.path.exists(file):
            path, fileN = os.path.split(file)
            os.system("(cd %s; zip -q %s.zip %s)" % (path, fileN, fileN))
            os.system("cp -u %s.zip %s" % (file, dir))
        else:
            print >>sys.stderr, "** %s does not exist" % file
#--------------------------------------------------------------


def copypng(dir, *fileL):
    for pngFile in fileL:
        if os.path.exists(pngFile):
            os.system("cp -u %s %s" % (pngFile, dir))
        else:
            print >>sys.stderr, "** %s does not exist" % pngFile
        if os.path.exists(pngFile.replace('png', 'pdf')):
            os.system("cp -u %s %s" % (pngFile.replace('png', 'pdf'), dir))
#--------------------------------------------------------------
def older(file, baseL):
    if not os.path.exists(file):
        return True
    file_mtime = time.localtime(os.stat(file).st_mtime)
    for base_file in baseL:
        base_file_mtime = time.localtime(os.stat(base_file).st_mtime)
        if file_mtime < base_file_mtime:
            return True
    return False
#----------------------------

def targetOlder(target, sourceL):
    if not os.path.exists(target):
        return True
    target_mtime = time.localtime(os.stat(target).st_mtime)
    for source_file in sourceL:
        source_file_mtime = time.localtime(os.stat(source_file).st_mtime)
        if target_mtime < source_file_mtime:
            return True
    return False
#----------------------------
def sourceNewer(source, targetL):
    source_mtime = time.localtime(os.stat(source).st_mtime)
    for target_file in targetL:
        if not os.path.exists(target_file):
            return True
        target_file_mtime = time.localtime(os.stat(target_file).st_mtime)
        if source_mtime > target_file_mtime:
            return True
    return False
#-----------------------------------------------------
#----------------------------

def targetNewer(file, baseL):
    return not targetOlder(file, baseL)


def copypdf(dir, *fileL):
    no_result = "/MPATHB/self/resource/sample/no_result.png"
    exist = 1
    for pdfFile in fileL:
        pngFile = pdfFile.replace('pdf', 'png')
        if os.path.exists(pdfFile):
            os.system("cp -u %s %s" % (pdfFile, dir))
            if not os.path.exists(pngFile) or \
                time.localtime(os.stat(pngFile).st_mtime) < \
                time.localtime(os.stat(pdfFile).st_mtime):
                os.system(' '.join([
                    "convert -density 150 -background white -alpha off",
                    pdfFile,pngFile]))
            os.system("cp -u %s %s" % (pngFile, dir))
        else:
            os.system("cp -u %s %s/%s" % (no_result, dir, os.path.split(pngFile)[1]))
            print >>sys.stderr, "Unexisted file %s" % pdfFile
            exist = 0
    return exist
#--------------------------------------------------------------
def knitr_read_txt(report_dir, curation_label):
    curation = 'CURATION_'+curation_label+'.txt'
    if not os.path.exists(report_dir+'/'+curation):
        fh = open(report_dir+'/'+curation, 'w')
        print >>fh, """
****************************

> CURATION
    
Add annotation here in markdown format for abnormal information if exists.
    
> CURATION

****************************    
"""
        fh.close()
    print """```{{r child="{}"}}
```
""".format(curation)
#-------------------------------------------
def getRelativeDir(file, report_sub_dir):
    if isinstance(file, list):
        tmpL = []
        for each_file in file:
            tmpL.append(report_sub_dir+'/'+os.path.split(each_file)[-1])
        return tmpL
    else:
        return report_sub_dir+'/'+os.path.split(file)[-1]
#--------------------------------------
def getFileName(file):
    if isinstance(file, list):
        tmpL = []
        for each_file in file:
            tmpL.append(os.path.split(each_file)[-1])
        return tmpL
    else:
        return os.path.split(file)[-1]
#--------------------------------------



def grenerateLabel():
    return time.strftime("L%Y-%m-%d-%H-%M-%S", time.localtime())



def generateLink(fileL, labelL, type, report_sub_dir='', join_symbol=' '):
    '''
    fileL: a list of files with relative path or absolute path
    labelL: a list of description names for each file in fileL
    type: normally pdf or xls to specify file type or extension
    report_sub_dir: a directory name; If given, relative path of fileL will be computed
    '''
    #Generate new relative location if `report_sub_dir` is given
    #Otherwise, we take files in fileL containing relative location
    if report_sub_dir:
        fileL = [report_sub_dir+'/'+os.path.split(file)[-1] \
                for file in fileL]
    #-------------------------------------
    #link = [label_type](file) [label2_type](file2)
    linkL = []
    for file, label in zip(fileL, labelL):
        tmp_71 = ['[', label, '_', type, ']', '(', file, ')']
        linkL.append(''.join(tmp_71))
    return join_symbol.join(linkL)

def grenerateQuotedLists(fileL, report_sub_dir='', quote_by="'"):
    '''
    fileL: a list of files with relative path or absolute path
    quote_by: single quote or double quotes
    report_sub_dir: a directory name; If given, relative path of fileL will be computed
    '''
    #Generate new relative location if `report_sub_dir` is given
    #Otherwise, we take files in fileL containing relative location
    if report_sub_dir:
        fileL = [report_sub_dir+'/'+os.path.split(file)[-1] \
                for file in fileL]
    #-------------------------------------
    #quotedLists = "'file1', 'file2', 'file3'"
    return ','.join([quote_by+file+quote_by for file in fileL])

def skipUnExistedOrEmptyFile(fileL, labelL=[]):
    newFileL = []
    newLabelL = []
    if not labelL:
        labelL = fileL
    for file, label in zip(fileL, labelL):
        if os.path.exists(file) and os.stat(file).st_size>0:
            newFileL.append(file)
            newLabelL.append(label)
    #----------------------------------------
    if labelL:
        return newFileL, newLabelL
    else:
        return newFileL


alphabet = [chr(i) for i in range(97,123)]

def num2alphabet(num):
    global alphabet
    division = num / 26
    remain = num % 26
    label = ''
    label = alphabet[remain] + label
    while division > 25:
        remain = division % 26
        division = division / 26
        label = alphabet[remain-1] + label
    #---------------------------
    if division > 0:
        label = alphabet[division-1] + label
    return label
#-----------------------------------


def checkLegalWord(word, stop=True):
    numstart = re.compile(r"^[0-9]")
    non_num_alphabet_underline = re.compile(r"[^0-9a-zA-Z_]")

    if numstart.search(word) or non_num_alphabet_underline.search(word):
        print >>sys.stderr, "Only numbers, alphabets, underlines are supported. And numbers should not be placed at first positions."
        if stop:
            sys.exit(1)
        return 1
    return 0
#-----checkLegalWord------------------


def transferListToMultiLTable(array, header=1, digit_format=1):
    '''
    This is designed to transfer 2-dimentional arrary to grid table
    for pandoc markdown.
    '''
    #---------Transfer digit ------------------
    maxLenL = [len(col) for col in array[0]]
    lenmaxLenL = len(maxLenL)
    alignL = ['%' if i.find('-') < 1 and \
            i.replace('.', '').replace('-', '').replace('e-','').isdigit() \
            else '%-' for i in array[1]]
    #print alignL
    for rowL in array[1:]:
        for i in range(lenmaxLenL):
            current = rowL[i]
            if current.find('-') < 1 and \
                current.replace('.', '').replace('-','').replace('e-','').isdigit() \
                and current != '0':
                if current.find('.') == -1:
                    current = "{:,}".format(int(current))
                else:
                    current = float(current)
                    if current < 1:
                        current = "{:,.4f}".format(current)
                    else:
                        current = "{:,.1f}".format(current)
                #------------------------------
                rowL[i] = current
            if len(current) > maxLenL[i]:
                maxLenL[i] = len(current)
    maxLenL = [i + 2 for i in maxLenL]
    #---------------------------------------- 
    gridL = []
    formatL = []
    for i in range(lenmaxLenL):
        maxLen = maxLenL[i]
        format = alignL[i] + str(maxLen) + 's'
        formatL.append(format)
    newarray = [[formatL[i] % rowL[i] for i in range(lenmaxLenL)] for
            rowL in array]
    if header:
        gridL.append("-"*70)
        gridL.append(' '.join(newarray[0]))
        gridL.append(' '.join(['-'*i for i in maxLenL]))
        newarray = newarray[1:]
    else:
        gridL.append(' '.join(['-'*i for i in maxLenL]))
    #------------------------
    for array88 in newarray[:-1]:
        gridL.append(' '.join(array88))
        gridL.append('')
    gridL.append(' '.join(newarray[-1]))
    if header:        
        gridL.append("-"*70)
    else:
        gridL.append(' '.join(['-'*i for i in maxLenL]))
    return gridL
#-----------------------------------------------
