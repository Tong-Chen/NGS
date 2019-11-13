#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Copyright 2018, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:

'''

import sys
import os
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import re
from collections import Counter
from itertools import islice
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print(desc, file=sys.stderr)
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="CDS fasta sequence.")
    parser.add_option("-o", "--output-prefix", dest="op",
        default="CRISPR", 
        metavar="FILEIN", help="Output prefix.")
    parser.add_option("-G", "--genome-index", dest="genome_index",
        default="/MPATHB/resource/Genome/DNA/Danshen/Danshen.bowtie", 
        metavar="FILEIN", help="Genome index for searching with bowtie. Default </MPATHB/resource/Genome/DNA/Danshen/Danshen.bowtie>.")
    parser.add_option("-b", "--bowtie-parameter", dest="bowtie_parameter",
        default="-v 0 -m 1", 
        metavar="FILEIN", help="Parameter for bowtie search. Default <-v 0 -m 1> meaning nomismatch only allow one match position in genome.")
    parser.add_option("-s", "--region-start", dest="region_start",
        type="int", default=200, help="Start site used for primer design. Default 200.")
    parser.add_option("-e", "--region-end", dest="region_end",
        type="int", default=700, help="End site used for primer design. Default 700.")
    parser.add_option("-p", "--primer-pattern", dest="primer_pattern",
        default='.{18}GG', help="Sequence pattern for primers. Default <.{18}GG> meaning primers should ends with GG and has length 20.")
    parser.add_option("-c", "--GC-percent", dest="GC_percent",
        default='50-70', help="GC content limitation for potential primers. Default <50-70>.")
    #parser.add_option("-l", "--primer-length", dest="primer_length",
    #    type="int", default=20, help="Primer length. Default 20.")
    parser.add_option("-g", "--guide-rna", dest="guide_rna",
        help="A fasta sequence specifying guide RNAs.")
    parser.add_option("-m", "--max-allowed-base-pair-with-guide-rna", dest="max_base_pair2guide",
        type="int", default=7, help="Max allowed base pair with guide RNA. Default 7. (<=7)")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

transD = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}

def readGuide(guide_file, transD):
    seqL = []
    aDict = {}
    for line in open(guide_file):
        if line[0] == '>':
            key =line.strip()[1:]
            aDict[key] = []
        else:
            aDict[key].append(line.strip())
    #------------------------------
    for key, valueL in aDict.items():
        seq = ''.join(valueL)
        seqL = [transD[i] for i in seq]
        seqL.reverse()
        seq = ''.join(seqL)
        aDict[key] = seq
    return aDict.values()
#---------------------------------

def generateKmers(seqL, k):
    kmer = set([])
    for seq in seqL:
        kmer.update(set([''.join(islice(seq, i, i+k)) for i in range(len(seq)+1-k)]))
    return kmer
#---------------------------------------------

def computeGC(seq, GC_percent):
    gc = int((seq.count('C') + seq.count('G')) * 100.0 / len(seq))
    if gc >= GC_percent[0] and gc <= GC_percent[1]:
        return gc
    else:
        return False

#---------------------------------------------------------
def match_guideRNA(seq, kmers, k):
    #for each_kmer in kmers:
    #    if seq.find(each_kmer) != -1:
    #        return True
    #return False
    seq_kmers = generateKmers([seq], k)
    if debug:
        print("\t\t\tKmers for potential primers", seq_kmers, file = sys.stderr)
    return seq_kmers.intersection(kmers)
    
#------------------------------------------------

def findPrimer(key, seq, kmers, max_base_pair2guide, primer_pattern, GC_percent, potential_primer_fh):
    if debug:
        print("\tIdentify primers for ", key, file=sys.stderr)
        print("\tIdentify primers for ", seq, file=sys.stderr)
    #seq_match = primer_pattern.finditer(seq)
    seq_match = primer_pattern.search(seq)
    primer_start = 0
    while seq_match:
        primer_start = primer_start + seq_match.span()[0]+1
        primer_seq   = seq_match.group()
        gc_ok = computeGC(primer_seq, GC_percent)
        if debug:
            print("\t\tPotential primers:", primer_seq, primer_start, file=sys.stderr)
        if gc_ok and (not match_guideRNA(primer_seq,  kmers, max_base_pair2guide)):
            seqname = ['>'+key, str(primer_start), str(gc_ok)]
            if debug:
                print("\t\tPotential primers OK:", primer_seq, file=sys.stderr)
            print("__".join(seqname), file=potential_primer_fh)
            print(primer_seq, file=potential_primer_fh)
        seq_match = primer_pattern.search(seq[primer_start:])
#------------------------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file2 = options.filein
    region_start = options.region_start - 1
    region_end   = options.region_end
    guide_rna    = options.guide_rna
    pattern      = options.primer_pattern
    primer_pattern = re.compile(pattern)
    GC_percent  = options.GC_percent
    GC_percent  = [float(i) for i in GC_percent.split('-')]
    #primer_length = options.primer_length
    max_base_pair2guide = options.max_base_pair2guide + 1
    op = options.op
    potential_primer = op + '.potential_primer.fa'
    potential_primer_fh = open(potential_primer, 'w')
    genome_index = options.genome_index
    bowtie_parameter = options.bowtie_parameter
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    kmers = set([])
    if guide_rna:
        guide_rnaL = readGuide(guide_rna, transD)
        kmers = generateKmers(guide_rnaL, max_base_pair2guide)
        if debug:
            print("Kmers for guide rna", guide_rnaL, file=sys.stderr)
            print(kmers, file=sys.stderr)

    if file2 == '-':
        fh = sys.stdin
    else:
        fh = open(file2)
    #--------------------------------
    count = 1
    key = ''
    for line in fh:
        if line[0] == '>':
            if key:
                seq = ''.join(seqL)[region_start:region_end]
                findPrimer(key, seq, kmers, max_base_pair2guide, primer_pattern, GC_percent, potential_primer_fh)
            key = line.split()[0][1:] + '__' + str(count)
            seqL = []
            count += 1
        else:
            seqL.append(line.strip())
    #-------------END reading file----------
    seq = ''.join(seqL)[region_start:region_end]
    findPrimer(key, seq, kmers, max_base_pair2guide, primer_pattern, GC_percent, potential_primer_fh)
    potential_primer_fh.close()
    #----close file handle for files-----
    if file2 != '-':
        fh.close()
    #-----------end close fh-----------

    if genome_index:
        newoutput = op+'.genome_one_match.fa'
        bowtie = ['bowtie -f --threads 30', bowtie_parameter, genome_index, potential_primer, "--al "+newoutput, '>/dev/null']
        bowtie = ' '.join(bowtie)
        if not os.system(bowtie):
            bowtie_parameter = newoutput

    final_output = op + '.crispr_primers.txt'
    final_output_fh = open(final_output, 'w')

    headerL = ["Gene", "Primer seq", "Primer start position", "Primer GC(%)", "Gene unique identifier (Genes order in input fasta file)"]
    print('\t'.join(headerL), file=final_output_fh)
    
    keyL = []
    for line in open(bowtie_parameter):
        if line[0] == '>':
            if keyL:
                tmpL = [keyL[0], seq, keyL[2], keyL[3], keyL[1]]
                print('\t'.join(tmpL), file=final_output_fh)
            keyL = line.strip()[1:].rsplit('__', 3)
        else:
            seq = line.strip()
    #----------------------------------------------------------
    if keyL:
        tmpL = [keyL[0], seq, keyL[2], keyL[3], keyL[1]]
        print('\t'.join(tmpL), file=final_output_fh)
    #----------------------------------------------------
    final_output_fh.close()
    
    
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
    if verbose:
        #print("--Successful %s" % strftime(timeformat, localtime()), file=sys.stderr)
        print >>sys.stderr, "--Successful %s" % strftime(timeformat, localtime())

if __name__ == '__main__':
    main()
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


