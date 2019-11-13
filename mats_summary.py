#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2013, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Functional description:
    1. Read the 'sumamry.txt' file from output of MATS and get the
       stacked bar plot.
    2. Select significant AS events.
    3. Get the sashimi plot of all significant AS events.
    3. Summary all AS events in one file, which can be given to
       <freeCompareWithName.py> to do comparison.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

debug = 0

def fprint(content):
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
    parser.add_option("-i", "--sample-file", dest="sampleFile",
        metavar="FILEIN", help="sampleFile given to mats.py")
    parser.add_option("-c", "--compare-pair", dest="filein",
        metavar="FILEIN", help="compare_pair given to mats.py")
    parser.add_option("-m", "--mats-result-dir", dest="mats_result_dir",
        help="The result dir given to mats.py ")
    parser.add_option("-g", "--genome-fa", dest="genome_fa",
        help="Genome fa file.")
    parser.add_option("-p", "--prefix", dest="prefix",
        default='mats', help="The prefix for output files.")
    parser.add_option("-q", "--fdr", dest="fdr",
        default=0.05, help="FDR for selecting significant AS")
    parser.add_option("-a", "--anno", dest="anno",
        help="Annotation file")
    parser.add_option("-s", "--sashimi", dest="sashimi",
        default=1, 
        help="Get the sashimi plot of significant AS events. \
Default <1> meaning get sashimi plots. \
Accept <0> to turn this off. \
Accept <2> to print sashimi plot command only.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readSummary(compare_pair, prefix, mats_result_dir):
    aDict = {}
    nameD = {'SE':'Exon skipping', 'MXE':'Mutually exclusive exon', 
            'RI':'Intron retention', 
            'A5SS':"Alternative 5' splice site", 
            'A3SS':"Alternative 3' splice site"}
    typeOrder = ','.join(['"' + "'Exon skipping'", 
        "'Intron retention'",
        "'Mutually exclusive exon'", 
        "'Alternative 5\\' splice site'", 
        "'Alternative 3\\' splice site'" +'"'])
    aDict['jcCount'] = {}
    aDict['sigJcCount'] = {}
    fileL = []
    for line in open(compare_pair):
        samp1, samp2 = line.split()
        file = '_vs_'.join([samp1, samp2])
        fileL.append(file)
        aDict['jcCount'][file] = {}
        aDict['sigJcCount'][file] = {}
        fh = open(mats_result_dir+'/'+file+'/summary.txt')
        line = fh.readline()
        while line.find("MATS Report") == -1:
            line = fh.readline()
        line = fh.readline() # ==== line
        line = fh.readline() # EventType NumEvents.JC.only  line
        line = fh.readline() # ==== line

        line = fh.readline()
        while line.find("===") == -1:
            lineL = line.split('\t')
            type  = lineL[0]
            jcCount = lineL[1]
            #print lineL[2]
            sigJcCount, sigJcCountDiff = lineL[2].split()
            sigJcUp, sigJcDw = sigJcCountDiff.strip(')(').split(':') 
            aDict['jcCount'][file][type] = jcCount    
            aDict['sigJcCount'][file][type] = {}
            aDict['sigJcCount'][file][type]['all'] = sigJcCount    
            aDict['sigJcCount'][file][type]['up'] = sigJcUp    
            aDict['sigJcCount'][file][type]['dw'] = sigJcDw    
            line = fh.readline()
    #-------END save--------------------------------------
    jcCount_f = 'MATS/' + prefix + '.AS.jcCount'
    sigJcCount_f = 'MATS/' + prefix + '.AS.sigJcCount'
    jcCount_fh = open(jcCount_f, 'w')
    sigJcCount_fh = open(sigJcCount_f, 'w')
    print >>jcCount_fh, "variable\tvalue\tSample"
    print >>sigJcCount_fh, "variable\tvalue\tSample"
    fileOrder = '"' + ','.join(['\''+file+'\'' for file in fileL]) + '"'
    for file in fileL:
        tmpDict = aDict['jcCount'][file]
        for type, value in tmpDict.items():
            print >>jcCount_fh, "\t".join([nameD[type], value, file]) 
    jcCount_fh.close()
    cmd = ['s-plot barPlot -f', jcCount_f, '-m TRUE -a Sample -l',
            typeOrder, "-L", fileOrder, '-R 90 -E pdf']
    if debug:
        print >>sys.stderr, ' '.join(cmd)
    os.system(' '.join(cmd))

    for file in fileL:
        tmpDict = aDict['sigJcCount'][file]
        for type, valueD in tmpDict.items():
            print >>sigJcCount_fh, "\t".join([nameD[type], valueD['all'], file]) 
    sigJcCount_fh.close()
    cmd = ['s-plot barPlot -f', sigJcCount_f, '-m TRUE -a Sample -l',
            typeOrder, "-L", fileOrder, '-R 90 -E pdf']
    if debug:
        print >>sys.stderr, ' '.join(cmd)
    os.system(' '.join(cmd))

    sigJcCount_f_ud = 'MATS/' + prefix + '.AS.sigJcCount.up_dw'
    sigJcCount_fh_ud = open(sigJcCount_f_ud, 'w')
    print >>sigJcCount_fh_ud, "variable\tvalue\tSample"
    fileL2 = []
    for file in fileL:
        tmpDict = aDict['sigJcCount'][file]
        for type, valueD in tmpDict.items():
            samp1, samp2 = file.split('_vs_')
            samp1 = file+'.'+samp1+'.up'
            samp2 = file+'.'+samp2+'.up'
            if samp1 not in fileL2:
                fileL2.append(samp1)
            if samp2 not in fileL2:
                fileL2.append(samp2)
            print >>sigJcCount_fh_ud, "\t".join([nameD[type],
                valueD['up'], samp1]) 
            print >>sigJcCount_fh_ud, "\t".join([nameD[type],
                valueD['dw'], samp2]) 
    sigJcCount_fh_ud.close()
    fileOrder = '"' + ','.join(['\''+file+'\'' for file in fileL2]) + '"'
    cmd = ['s-plot barPlot -f', sigJcCount_f_ud, '-m TRUE -a Sample -l',
            typeOrder, "-L", fileOrder, '-R 90 -E pdf']
    if debug:
        print >>sys.stderr, ' '.join(cmd)
    os.system(' '.join(cmd))
#-----------------------------------
def getReverseComplement(seq, type='DNA'):
    if type == 'DNA':
        dict = {'A':'T', 'G':'C', 'T':'A','C':'G','a':'t','g':'c','t':'a','c':'g', 'N':'N', 'n':'n'}
    elif type == 'RNA':
        dict = {'A':'U', 'G':'C', 'U':'A','C':'G','a':'u','g':'c','u':'a','c':'g'}
    #-------------------------------------------
    seqL = list(seq)
    seqL.reverse()
    return ''.join([dict[i] for i in seqL])
#----------------------------------------------------

def getseq(lineL, type, genomeD):
    if not genomeD:
        return '-', '-'
    chr = lineL[3]
    chrseq = genomeD[chr]
    strand = lineL[4]
    if type in ["A3SS", "A5SS"]:
        flankingES = int(lineL[9])
        flankingEE = int(lineL[10])
        common = chrseq[flankingES:flankingEE]
        if strand == '-':
            common = getReverseComplement(common)
        common = common.lower()
        longExonStart_0base = int(lineL[5])
        longExonEnd = int(lineL[6])
        AsSeq = chrseq[longExonStart_0base:longExonEnd]
        if strand == '-':
            AsSeq = getReverseComplement(AsSeq)
        shortES = int(lineL[7])
        shortEE = int(lineL[8])
        NormalSeq = chrseq[shortES:shortEE]
        if strand == '-':
            NormalSeq = getReverseComplement(NormalSeq)
        if type == "A3SS":
            AsSeq, NormalSeq = common+AsSeq, common+NormalSeq
        elif type == "A5SS":
            AsSeq, NormalSeq = AsSeq+common, NormalSeq+common
    elif type in ["SE", "RI"]:
        upstreamES = int(lineL[7])
        upstreamEE = int(lineL[8])
        up_t = chrseq[upstreamES:upstreamEE]
        if strand == '-':
            dw = getReverseComplement(up_t)
            dw = dw.lower()
        else:
            up = up_t.lower()
        downstreamES = int(lineL[9])
        downstreamEE = int(lineL[10])
        dw_t = chrseq[downstreamES:downstreamEE]
        if strand == '-':
            up = getReverseComplement(dw_t)
            up = up.lower()
        else:
            dw = dw_t.lower()
        exonStart_0base = int(lineL[5])
        exonEnd = int(lineL[6])
        if type == "SE":
            se = chrseq[exonStart_0base:exonEnd]
            if strand == '-':
                se = getReverseComplement(se)
            AsSeq, NormalSeq = up+se+dw, up+dw
        elif type == "RI":
            ri = chrseq[exonStart_0base:exonEnd]
            if strand == '-':
                ri = getReverseComplement(ri)
            ri = ri.replace(up.upper(), up)
            ri = ri.replace(dw.upper(), dw)
            AsSeq, NormalSeq = ri, up+dw
    elif type == "MXE":
        stExonStart_0base = int(lineL[5])
        stExonEnd         = int(lineL[6])
        st = chrseq[stExonStart_0base:stExonEnd]
        if strand == '-':
            st = getReverseComplement(st)
        ndExonStart_0base = int(lineL[7])
        ndExonEnd         = int(lineL[8])
        nd = chrseq[ndExonStart_0base:ndExonEnd]
        if strand == '-':
            nd = getReverseComplement(nd)
        upstreamES        = int(lineL[9])
        upstreamEE        = int(lineL[10])
        up_t = chrseq[upstreamES:upstreamEE]
        if strand == '-':
            dw = getReverseComplement(up_t)
            dw = dw.lower()
        else:
            up = up_t.lower()
        downstreamES      = int(lineL[11])
        downstreamEE      = int(lineL[12])
        dw_t = chrseq[downstreamES:downstreamEE]
        if strand == '-':
            up = getReverseComplement(dw_t)
            up = up.lower()
        else:
            dw = dw_t.lower()
        AsSeq,  NormalSeq = up+st+dw, up+nd+dw
        #-------------------------------------------
    #-------------------------------------------
    return AsSeq, NormalSeq
#-----------------------------------


def selectSig(compare_pair, prefix, fdr, sampD, sashimi, mats_result_dir, annoD={}, headerline='', genomeD={}):
    #print annoD['VIT_13s0084g00310']
    #print annoD['VIT_12s0034g02350']
    typeL = ['SE', 'MXE', 'A5SS', 'A3SS', 'RI']
    suffix = '.MATS.JunctionCountOnly.txt'
    #fileL = [i+suffix for i in fileL]
    sigFileD = {}

    for line in open(compare_pair):
        samp1, samp2 = line.strip().split('\t')
        folder = '_vs_'.join([samp1, samp2])
        rep1_bam = ','.join([rep+'/'+rep+'.Aligned.sortedByCoord.out.bam' for rep in
            sampD[samp1]])
        rep2_bam = ','.join([rep+'/'+rep+'.Aligned.sortedByCoord.out.bam' for rep in
            sampD[samp2]])
        for type in typeL:
            file = type + suffix
            anno_file = 'MATS/' + folder + '.' + file + '.anno.xls'
            anno_file = open(anno_file, 'w')
            samp1_up = folder + '.' + samp1 + '.UP'
            samp1_sig_anno = 'MATS/' + samp1_up + '.' + \
                str(fdr)+'.' + file + '.anno.xls'
            samp1_sig_anno = open(samp1_sig_anno, 'w')
            samp2_up = folder + '.' + samp2 + '.UP'
            samp2_sig_anno = 'MATS/' + samp2_up + '.' + \
                str(fdr)+'.'+ file + '.anno.xls'
            samp2_sig_anno = open(samp2_sig_anno, 'w')

            samp1_sig = 'MATS/' + samp1_up + '.' + \
                str(fdr)+'.' + file
            if samp1_up not in sigFileD:
                sigFileD[samp1_up] = []
            sigFileD[samp1_up].append([samp1_sig, type])
            samp1_sig_ = open(samp1_sig, 'w')
            samp2_sig = 'MATS/' + samp2_up + '.' + \
                str(fdr)+'.'+ file
            if samp2_up not in sigFileD:
                sigFileD[samp2_up] = []
            sigFileD[samp2_up].append([samp2_sig, type])

            samp2_sig_ = open(samp2_sig, 'w')

            file = mats_result_dir + '/' + folder + '/MATS_output/' + file
            header = 1
            for line in open(file):
                line = line.strip()
                if header:
                    as_header = line + "\tAsSeq\tNormalSeq"
                    as_header = as_header.replace('SAMPLE_1', samp1)
                    as_header = as_header.replace('SAMPLE_2', samp2)
                    if annoD:
                        print >>anno_file, "%s\t%s" % (as_header,headerline)
                        print >>samp1_sig_anno, "%s\t%s" % (as_header,headerline)
                        print >>samp2_sig_anno, "%s\t%s" % (as_header,headerline)
                        print >>samp1_sig_, as_header
                        print >>samp2_sig_, as_header
                    else:
                        print >>anno_file, as_header
                        print >>samp1_sig_anno, as_header
                        print >>samp2_sig_anno, as_header
                        print >>samp1_sig_, as_header
                        print >>samp2_sig_, as_header
                    #---------------------------------
                    header -= 1
                    continue
                #--------------------------------
                lineL = line.split('\t')
                gene = lineL[1].split('__')[0].replace("\"", '')
                anno = annoD.get(gene, "")
                fdr_cur  = float(lineL[19])
                if type == 'MXE':
                    fdr_cur = float(lineL[21])
                incLevelD = float(lineL[-1])
                AsSeq, NormalSeq = getseq(lineL, type, genomeD)
                print >>anno_file, "%s\t%s\t%s\t%s" % (line, AsSeq, NormalSeq, anno)
                if fdr_cur <= fdr:
                    if incLevelD > 0: #samp1
                        print >>samp1_sig_anno, "%s\t%s\t%s\t%s" % (line, AsSeq, NormalSeq, anno)
                        print >>samp1_sig_, line
                    elif incLevelD < 0: #samp2
                        print >>samp2_sig_anno, "%s\t%s\t%s\t%s" % (line, AsSeq, NormalSeq, anno)
                        print >>samp2_sig_, line
                #--------------------------------------
            #----------------------------------------
            anno_file.close()
            samp1_sig_anno.close()
            samp2_sig_anno.close()
            samp1_sig_.close()
            samp2_sig_.close()
            if sashimi:
                cmd = ["rmats2sashimiplot -b1", rep1_bam, "-b2",
                        rep2_bam, "-t", type, "-e", samp1_sig, "-l1",
                        samp1, "-l2", samp2, "-exon_s 1 -intron_s 1",
                        "-o", samp1_sig.replace('txt', "sashimi"),
                        '>', samp1_sig.replace('txt', "sashimi.log"), 
                        '2>&1 &']
                if sashimi == 2:
                    print " ".join(cmd)
                else:
                    os.system(" ".join(cmd))
                cmd = ["rmats2sashimiplot -b1", rep1_bam, "-b2",
                        rep2_bam, "-t", type, "-e", samp2_sig, "-l1",
                        samp1, "-l2", samp2, "-exon_s 1 -intron_s 1",
                        "-o", samp2_sig.replace('txt', "sashimi"),
                        '>', samp2_sig.replace('txt', "sashimi.log"), 
                        '2>&1 &']
                if sashimi == 2:
                    print " ".join(cmd)
                else:
                    os.system(" ".join(cmd))
        #---END one type---------------------------------------------
    #------------------------------------------------------           
    return sigFileD
#---------------selectSig------------------------------

def readMATS_output(file, event_type):
    asD = {}
    for line in open(file):
        geneSymbol = "";
        gene_no_str = "";
        id_str = "";
        if line.startswith('ID'):
          continue;
        items = line.split("\t")
        geneSymbol = items[2]
        geneSymbol = geneSymbol.replace('\"','');
        chr = items[3]
        strand = items[4]
        e1st_s = "";
        e1st_e = "";
        e2st_s = "";
        e2st_e = "";
        se_s = "";
        se_e = "";
        up_s = "";
        up_e = "";
        dn_s = "";
        dn_e = "";
        inc_level1 = "";
        inc_level2 = "";
        if event_type!='MXE':
            se_s = str(int(items[5])+1)
            se_e = items[6]
            up_s = str(int(items[7])+1)
            up_e = items[8]
            dn_s = str(int(items[9])+1)
            dn_e = items[10]
            inc_level1 = items[20] ## IncLevel1
            inc_level2 = items[21] ## IncLevel2
        if event_type=='MXE':
            e1st_s = str(int(items[5])+1)
            e1st_e = items[6]
            e2st_s = str(int(items[7])+1)
            e2st_e = items[8]
            up_s = str(int(items[9])+1)
            up_e = items[10]
            dn_s = str(int(items[11])+1)
            dn_e = items[12]
        #--------------------------------------------
        if strand =='+':
            #name_str = geneSymbol+"___"+chr+"_"+up_s+"_"+up_e+"_"+strand+"@"+chr+"_"+se_s+"_"+se_e+"_"+strand+"@"+chr+"_"+dn_s+"_"+dn_e+"_"+strand
            if event_type != 'MXE':
                id_str = "_".join([chr, up_s, up_e,
                    strand+"@"+chr, se_s, se_e, strand+"@"+chr,
                    dn_s, dn_e, strand])
            else:
                id_str = "_".join([chr, up_s, up_e,
                    strand+"@"+chr, e1st_s, e1st_e, 
                    strand+"@"+chr, e2st_s, e2st_e, 
                    strand+"@"+chr, dn_s, dn_e, strand])
            name_str = geneSymbol+"___"+id_str
        elif strand =='-':
            if event_type != 'MXE':
                id_str = "_".join([chr, dn_s, dn_e,
                    strand+"@"+chr, se_s, se_e, 
                    strand+"@"+chr, up_s, up_e, strand])
            else:
                id_str = "_".join([chr, dn_s, dn_e,
                    strand+"@"+chr, e2st_s, e2st_e, 
                    strand+"@"+chr, e1st_s, e1st_e, 
                    strand+"@"+chr, up_s, up_e, strand])
            name_str = geneSymbol+"___"+id_str
        assert name_str not in asD, "Duplictae "+ name_str + '...' + file
        asD[name_str] = line
    #---------------------END
    return asD
#-----------------------------------------

def getTotal(sigFileD, prefix):
    '''
    sigFileD = {  
            samp1_up: [[samp1_sig, 'SE'], [samp1_sig, 'RI']], 
            samp2_up: [[samp2_sig, 'SE'], [samp2_sig, 'RI']]
        } 
    '''
    outputAll = open('MATS/'+prefix+'.allType.DE', 'w')
    typeL = ['SE', 'MXE', 'A5SS', 'A3SS', 'RI']
    outputFD = dict([[type, 'MATS/'+prefix+'.'+type+'.DE'] for type in typeL]) 
    #print outputFD
    outputFH = dict([[type, open(outputFD[type], 'w')] for type in typeL])
    #print sigFileD
    for compare, typeL in sigFileD.items():
        for file, type in typeL:
            #print type
            asD = readMATS_output(file, type)
            for key in asD.keys():
                #print "%s\t%s" % (key, compare)
                print >>outputFH[type], "%s\t%s" % (key, compare)
                gene, pos = key.split("___")
                print >>outputAll, "\t".join([gene, pos, type, compare])
        #----------END each type-----------------
    #----------------------------------------------
    outputAll.close()
    for fileH in outputFH.values():
        fileH.close()
#------getTotal------------------

def readAnno(anno):
    header = 1
    annoD = {}
    for line in open(anno):
        if header:
            headerline = line.strip()
            header -= 1
            continue
        #----------------------------
        line = line.strip()
        lineL = line.split('\t', 1)
        key = lineL[0]
        annoD[key] = line
    #-----------------------------------
    #print annoD.keys()[1:10]
    #print annoD['VIT_13s0084g00310']
    #print annoD['VIT_12s0034g02350']
    return annoD, headerline
#---------------------------------------
def readGenome(genome_fa):
    genomeD = {}
    for line in open(genome_fa):
        if line[0] == '>':
            key = line.split()[0][1:]
            assert key not in genomeD, "Duplicate key "+key
            genomeD[key] = []
        else:
            genomeD[key].append(line.strip().upper())
    #-------------------------------------------
    for key, seqL in genomeD.items():
        genomeD[key] = ''.join(seqL)
    return genomeD
#------------readGenome-----------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    compare_pair = options.filein
    sampleFile = options.sampleFile
    mats_result_dir = options.mats_result_dir
    sampleD = {}
    for line in open(sampleFile):
        rep, samp = line.strip().split()[:2]
        if samp not in sampleD:
            sampleD[samp] = []
        sampleD[samp].append(rep)
    #----------------------------
    genome_fa = options.genome_fa
    if genome_fa:
        genomeD = readGenome(genome_fa)
    else:
        genomeD = {}
    prefix = options.prefix
    fdr = float(options.fdr)
    anno = options.anno
    sashimi = int(options.sashimi)
    annoD = {}
    headerline = ''
    if anno:
        annoD, headerline = readAnno(anno)
    #--------------------------------------
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    readSummary(compare_pair, prefix, mats_result_dir)
    sigFileD = selectSig(compare_pair, prefix, fdr, sampleD, sashimi, mats_result_dir, annoD, headerline, genomeD)
    #print sigFileD
    getTotal(sigFileD, prefix)
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


