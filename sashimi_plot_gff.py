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

Input GFF for one gene or whole GFF file. The program will selcect exons within given coordinates.

1. Only <exon> in third column will be used.
2. Different <transcript_id> in the ninth column will be used to represent different transcripts. All other fields will be ignored.

4	Cufflinks	exon	921517	922929	.	-	.	gene_id "XLOC_021911"; transcript_id "grape_00072032:VIT_04s0008g01110:HSFA2"; exon_number "1"; gene_name "gene:VIT_04s0008g01110.1"; oId "T45_3.16011.2"; nearest_ref "transcript:VIT_04s0008g01110.t01"; class_code "j"; tss_id "TSS27177";
4	Cufflinks	exon	923502	923881	.	-	.	gene_id "XLOC_021911"; transcript_id "grape_00072032:VIT_04s0008g01110:HSFA2"; exon_number "2"; gene_name "gene:VIT_04s0008g01110.1"; oId "T45_3.16011.2"; nearest_ref "transcript:VIT_04s0008g01110.t01"; class_code "j"; tss_id "TSS27177";
4	Cufflinks	exon	921517	921930	.	-	.	gene_id "XLOC_021911"; transcript_id "grape_00072031:VIT_04s0008g01110:HSFA2"; exon_number "1"; gene_name "gene:VIT_04s0008g01110.1"; oId "T45_3.16011.1"; nearest_ref "transcript:VIT_04s0008g01110.t01"; class_code "o"; tss_id "TSS27177";
4	Cufflinks	exon	922020	923881	.	-	.	gene_id "XLOC_021911"; transcript_id "grape_00072031:VIT_04s0008g01110:HSFA2"; exon_number "2"; gene_name "gene:VIT_04s0008g01110.1"; oId "T45_3.16011.1"; nearest_ref "transcript:VIT_04s0008g01110.t01"; class_code "o"; tss_id "TSS27177";
4	Cufflinks	exon	921539	921930	.	-	.	gene_id "XLOC_021911"; transcript_id "grape_00072755:VIT_04s0008g01110:HSFA2"; exon_number "1"; gene_name "gene:VIT_04s0008g01110.1"; oId "T40_2.16151.2"; nearest_ref "transcript:VIT_04s0008g01110.t01"; class_code "j"; tss_id "TSS27178";
4	Cufflinks	exon	922020	922929	.	-	.	gene_id "XLOC_021911"; transcript_id "grape_00072755:VIT_04s0008g01110:HSFA2"; exon_number "2"; gene_name "gene:VIT_04s0008g01110.1"; oId "T40_2.16151.2"; nearest_ref "transcript:VIT_04s0008g01110.t01"; class_code "j"; tss_id "TSS27178";
4	Cufflinks	exon	923502	923728	.	-	.	gene_id "XLOC_021911"; transcript_id "grape_00072755:VIT_04s0008g01110:HSFA2"; exon_number "3"; gene_name "gene:VIT_04s0008g01110.1"; oId "T40_2.16151.2"; nearest_ref "transcript:VIT_04s0008g01110.t01"; class_code "j"; tss_id "TSS27178";
4	Cufflinks	exon	948963	948983	.	-	.	gene_id "XLOC_021911"; transcript_id "grape_00072755:VIT_04s0008g01110:HSFA2"; exon_number "4"; gene_name "gene:VIT_04s0008g01110.1"; oId "T40_2.16151.2"; nearest_ref "transcript:VIT_04s0008g01110.t01"; class_code "j"; tss_id "TSS27178";

Output GFF will be saved in file <name_sashimiplot/name.gff3> (name is given in coordinates)

* One gene line with ID and Name
* Several mRNA lines with individual ID and common gene parent
* Several exon lines for each mRNA with individual ID and common mRNA parent

4	ct	gene	921517	923881	.	-	.	ID=HSFA2;Name=HSFA2
4	ct	mRNA	921517	923881	.	-	.	ID=grape_00072031:VIT_04s0008g01110:HSFA2;Parent=HSFA2
4	ct	exon	921517	921930	.	-	.	ID=1;Parent=grape_00072031:VIT_04s0008g01110:HSFA2
4	ct	exon	922020	923881	.	-	.	ID=2;Parent=grape_00072031:VIT_04s0008g01110:HSFA2
4	ct	mRNA	921539	923728	.	-	.	ID=grape_00072755:VIT_04s0008g01110:HSFA2;Parent=HSFA2
4	ct	exon	921539	921930	.	-	.	ID=1;Parent=grape_00072755:VIT_04s0008g01110:HSFA2
4	ct	exon	922020	922929	.	-	.	ID=2;Parent=grape_00072755:VIT_04s0008g01110:HSFA2
4	ct	exon	923502	923728	.	-	.	ID=3;Parent=grape_00072755:VIT_04s0008g01110:HSFA2

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import re
comma_semicolon = re.compile(r'[;,]')
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
        metavar="FILEIN", help="A gff3 or gtf file with format specified above. Multiple genes can be included.")
    parser.add_option("-c", "--coordinate", dest="coordinate",
            help="Supply gene coordinate in format like <chr:start:end:strand:name>. \
Multiple coordinates separated by ';' are accepted. [The meaning of <start> and <end> ] \
is same as in gff file. 1-based, end not included.")
    parser.add_option("-b", "--bamFiles", dest="bamFiles",
        help="BAM files to be used in format like <samp1_1.bam,samp1_2.bam;samp2_1.bam;samp3_1.bam,samp3_2.bam,samp3_3.bam>. ")
    parser.add_option("-l", "--bamLabels", dest="bamLabels",
        help="Labels for BAM files given to -b in format like <samp1_1,samp1_2;samp2_1;samp3_1,samp3_2,samp3_3>. ")
    parser.add_option("-r", "--total-mapped-reads-count", dest="readsCnt",
        help="Total mapped reads count used for normalize BAM counts in format like <10000,10001;10010;10020,10021,10019>. ")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def read(gff, coordinate):
    '''
    gffD = {tr_id1: [
                [chr, ct, exon, start, end, ., strand, .], 
                [chr, ct, exon, start, end, ., strand, .]], 
            tr_id2: [
                [chr, ct, exon, start, end, ., strand, .]]
            }

    mRNAD = {tr_id1: '\t'.join([chr, ct, mRNA, start, end, ., strand, .])}
    '''
    chr, start, end, strand, name = coordinate.split(':')
    start = int(start)
    end   = int(end)
    gffD = {}
    for line in open(gff):
        if line[0] == '#':
            continue
        lineL = line.strip().strip(';').split('\t')
        chr_this = lineL[0]
        lineL[1] = 'ct'
        type = lineL[2]
        if chr_this == chr and type == 'exon':
            exon_s = int(lineL[3])
            exon_e = int(lineL[4])
            exon_strand = lineL[6]
            if exon_strand == strand and (start<=exon_s) and (end>=exon_e):
                if debug:
                    print >>sys.stderr, lineL[8]
                annoD = dict([i.strip().split(' ') for i in lineL[8].replace('"', '').split(';')])
                tr_id = annoD['transcript_id']
                if tr_id not in gffD:
                    gffD[tr_id] = [lineL[:8]]
                else:
                    gffD[tr_id].append(lineL[:8])
        #---------------------------------------------
    #---------------------------------
    mRNAD = {}
    for tr_id, posL in gffD.items():
        posL.sort(key=lambda x: int(x[3]))
        start = posL[0][3]
        end   = posL[-1][4]
        mRNAD[tr_id] = '\t'.join([chr, 'ct\tmRNA', start, end, '.', strand, '.'])
    return gffD, mRNAD           
            

#----------------------------------
def transferGFF(gff, coordinate):
    '''
    coordinate: specify gene positon in format like <chr:start:end:strand:name>
    '''
    chr, start, end, strand, name = coordinate.split(':')
    aDict = {'chr':chr, 'start':start, 'end':end, 'strand':strand, 'id':name, 'name':name}
    dir = name+'_sashimiplot/'
    os.system("mkdir -p "+dir)
    gff_output = dir + name + '.gff3'
    fh = open(gff_output, 'w')
    print >>fh, "{d[chr]}\tct\tgene\t{d[start]}\t{d[end]}\t.\t{d[strand]}\t.\tID={d[id]};Name={d[name]}".format(d=aDict)
    #print >>fh, "{d[chr]}\tct\tmRNA\t{d[start]}\t{d[end]}\t.\t{d[strand]}\t.\tID=ct0;Name=ct0;Parent={d[name]}".format(d=aDict)
    gffD, mRNAD = read(gff, coordinate)
    for tr_id, mRNA in mRNAD.items():
        print >>fh, '{}\tID={};Parent={}'.format(mRNA, tr_id, name)
        count = 1
        for exonL in gffD[tr_id]:
            print >>fh, '{}\tID={};Parent={}'.format('\t'.join(exonL), count, tr_id)
            count += 1
        #------------------------------------
    fh.close()
    return gff_output, name, dir
#-------------------------------
def hex2rgb(hexcolor):
    return [(hexcolor>>16) & 0xff, (hexcolor>>8) & 0xff, hexcolor & 0xff]

def rgb2hex(rgbcolor):
    r, g, b = rgbcolor
    rgb = hex((r << 16) + (g << 8) +b)[2:].upper()
    zero = '0'* (6-len(rgb))
    return '#'+zero+rgb
#----------------------------------
def generateColor(labelL):
    colorL = []
    r = 255
    g = 255
    b = 255
    len_label = int(len(labelL) / 3 + 1)
    step = int(250 / len_label)
    
    cnt = 1
    for labels in labelL:
        if cnt % 3 == 1:
            r = r - step 
        elif cnt % 3 == 2:
            g = g -step
        else:
            b = b - step
        cnt += 1
        color = rgb2hex([r, g, b])
        for label in labels.split(','):
            colorL.append('"'+color+'"')
    return ','.join(colorL)

def generateSeeting(bamL, labelL, readsCnt, name, dir):
    #bam = ','.join(['"'+bams+'"' for bamS in bamL for bams in bamS.split(',')])
    bam = ','.join(['"'+bams+'"' for bams in comma_semicolon.split(bamL)])
    bamCnt = bam.count(',')+1
    #label = ','.join(['"'+labels+'"' for labelS in labelL for labels in labelS.split(',')])
    label = ','.join(['"'+labels+'"' for labels in comma_semicolon.split(labelL)])
    labelCnt = label.count(',')+1
    assert bamCnt == labelCnt, "{}\n{}".format(bamL, labelL)
    if readsCnt:
        coverage = ','.join(comma_semicolon.split(readsCnt))
    else:
        coverage = ','.join([1000000000]*labelCnt)
    coverageCnt = coverage.count(',')+1
    assert bamCnt == coverageCnt, "{}\n{}".format(bamL, readsCnt)
    
    color = generateColor(labelL.split(';'))
    setting = dir+'sashimi_plot_settings.txt'
    set_fh = open(setting, 'w')
    print >>set_fh, '''[data]
bam_prefix = 
miso_prefix = 
bam_files = [{d[bam]}]
miso_files = [{d[bam]}]
[plotting]
fig_width = 10
fig_height = 7
exon_scale = 1
intron_scale = 1
logged = False
font_size = 7
bar_posteriors = False
nyticks = 4
nxticks = 11
show_ylabel = True
show_xlabel = True
plot_title = "{d[name]}"
plot_label = plot_label
show_posteriors = False
number_junctions = True
resolution = .5
colors = [{d[color]}]
coverage = [{d[coverage]}]
sample_labels = [{d[label]}]
reverse_minus = True
    '''.format(d={"bam":bam, "label":label, 'name':name, 'color':color, 'coverage':coverage})
    set_fh.close()
    return setting
#-----------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    gff = options.filein
    coordinateL = options.coordinate.split(';')
    bamL = options.bamFiles
    labelL = options.bamLabels
    readsCnt = options.readsCnt
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    for coordinate in coordinateL:
        gff_output, name, dir = transferGFF(gff, coordinate)
        setting = generateSeeting(bamL, labelL, readsCnt, name, dir)
        index_path = dir+'sashimi_gff_index/'
        if os.path.exists(index_path):
            os.system("/bin/rm -rf "+index_path)
        cmd = ["index_gff --index", gff_output, index_path]
        
        if debug:
            print >>sys.stderr, " ".join(cmd)
        else:
            os.system(" ".join(cmd))
        cmd = ['sashimi_plot --plot-event', name, index_path, setting, "--output-dir", dir ]
        if debug:
            print >>sys.stderr, " ".join(cmd)
        else:
            os.system(" ".join(cmd))
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


