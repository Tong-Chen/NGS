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
    This is first designed to perform multiple sequence alignment and
    trimming for protein sequence and coding sequences using `mafft`
    and `trimal`.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
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
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-dir", dest="dir",
        help="A directory containing multiple \
pairs of protein sequences and coding sequences. All files \
are in Fasta format. The name of protein sequences and coding \
sequences for each pair. Protein file should end with `pep.fa` \
and cds file ends with `nucl.fa`. Normally the output of \
`parseOrthoMclResult.py`.")
    parser.add_option("-m", "--mafft-parameter", dest="mafftP",
        default="--maxiterate 1000 --genafpair --thread 30", 
        help="Parameters given to multiple sequence alignment \
software `mafft`. Default <--maxiterate 1000 --genafpair --thread 30>.")
    parser.add_option("-t", "--trimal-parameter", dest="trimalP",
        default="-automated1", 
        help="Parameters given to multiple sequence alignment trimming \
software `trimal`. Default <-automated1>.")
    parser.add_option("-C", "--COG", dest="cog",
        default='auto_CT', 
        help="COG file used fot `ete build`. For orthoMCL output \
The program will generate COG file automatically. For \
other input files, a COG file should be supplied.")
    parser.add_option("-c", "--clearall", dest="clearall",
        default=False, action="store_true", 
        help="Clear existing results.")
    parser.add_option("-p", "--phylogenetic-tree", dest="ptree",
        default='raxml_default_bootstrap', 
        help="Method to construct phylogenetic tree, \
accept <fasttree_default>, <raxml_default_bootstrap(default)> \
<phyml_default_bootstrap>, <phyml_default>, <raxml_default>.")
    parser.add_option("-o", "--output-prefix", dest="output_prefix",
        help="A string to label output files")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.dir != None, "A dir-name needed for -i"
    return (options, args)
#--------------------------------------------------------------------
def concatenate(trim_pepL, trim_pep_concat_file):
    '''
    file content:
    >C18780_And12391.1
    MAGTSLMDSLFQRSLDDLIKGIRLCPPGTEPAFI....
    >C18780_Sly03g119270.1.1
    MAGPSLLDSLFQRSLEDLIKGLRLFVG
    '''
    fh = open(trim_pep_concat_file, 'w')
    for trim_pep_file in trim_pepL:
        for line in open(trim_pep_file):
            if line[0] == '>':
                seqname, spename = line[1:].split('|', 1)[0].split('_')
                print >>fh, '>%s_%s' % (spename, seqname)
            else:
                print >>fh, line,
    #------------------------------------------
    fh.close()

#------------concatenate----------------------

def generateCOG(fileL, cog_file):
    '''
    file content:
    >C18780_And12391.1
    MAGTSLMDSLFQRSLDDLIKGIRLCPPGTEPAFI....
    >C18780_Sly03g119270.1.1
    MAGPSLLDSLFQRSLEDLIKGLRLFVG

    '''
    fh = open(cog_file, 'w')
    for file in fileL:
        spenameL = []
        for line in open(file+'.pep.fa'):
            if line[0] == '>':
                seqname, spename = line[1:].split('|', 1)[0].split('_')
                spenameL.append(spename)
        #--------------------------------
        spenameL = [spe+'_'+seqname for spe in spenameL]
        spenameL.sort()
        print >>fh, '\t'.join(spenameL)
    fh.close()
#----------generateCOG--------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    dir = options.dir
    mafftP = options.mafftP
    trimalP = options.trimalP
    cog_file = options.cog
    clearall = options.clearall
    ptree = options.ptree
    output_prefix = options.output_prefix+'.'+ptree
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    
    fileL = [dir+'/'+file.replace('.pep.fa', '') \
        for file in os.listdir(dir) if file.endswith('pep.fa')]

    if cog_file == 'auto_CT':
        cog_file = output_prefix + '.cog'
        generateCOG(fileL, cog_file)

    trim_pepL = []
    trim_nuclL = []
    for file in fileL:
        out_pep_aln = file+'.pep_mafft_einsi.aln'
        msa = ' '.join(['mafft', mafftP, file+'.pep.fa',
                        '>'+out_pep_aln, "2>"+file+'.log'])

        trim_pep_file = file+'.prot_mafft_einsi.trim.aln'
        trim_pepL.append(trim_pep_file)
        #trim_pep = ' '.join(['trimal -phylip -in', out_pep_aln, '-out',
        #        trim_pep_file, trimalP, "2>>"+file+'.log'])
        trim_pep = ' '.join(['trimal -in', out_pep_aln, '-out',
                trim_pep_file, trimalP, "2>>"+file+'.log'])
        if debug:
            print >>sys.stderr, msa
            print >>sys.stderr, trim_pep
        else:
            if clearall or not os.path.exists(out_pep_aln):
                os.system(msa)
            if verbose:
                print >>sys.stderr, "**MSA finished for %s" % file
            if clearall or not os.path.exists(trim_pep_file):
                os.system(trim_pep)
            if verbose:
                print >>sys.stderr, "**MSA trimming finished for %s" % file
        #-------------------------------------
        #----------backtrans--------------------
        if os.path.exists(file+'.nucl.fa'):
            trim_nucl_file = file+'.nucl_mafft_einsi.trim.aln'
            trim_nuclL.append(trim_nucl_file)
            trim_nucl = ' '.join(['trimal -in', out_pep_aln, '-out',
                    trim_nucl_file, trimalP, 
                    "-backtrans", file+'.nucl.fa', '2>>'+file+'.log'])
            if debug:
                print >>sys.stderr, trim_nucl
            else:
                if clearall or not os.path.exists(trim_nucl_file):
                    os.system(trim_nucl)
                if verbose:
                    print >>sys.stderr, "**Protein MSA based CDS alignment and trimming finished for %s\n" % file
        
        #-----------------------------------------
    #------------END align-trim-backtrans all files----------
    #----Concatenate trimmed files-------
    trim_pep_concat_file = output_prefix + '.trim_concat_prot.fa'

    concatenate(trim_pepL, trim_pep_concat_file)
            
    tree_build = ' '.join(["xvfb-run ete3 build --cpu 30 -w none-none-none-none", 
            "-m cog_100-alg_concat_default-%s" % ptree, 
            '-o %s' % output_prefix+'.sptree', '--clearall', 
            '-a %s' % trim_pep_concat_file, '--cogs %s' % cog_file, 
            '--tools-dir /root/.etetoolkit/ext_apps-latest/'])
    if debug:
        print >>sys.stderr, tree_build
    else:
        os.system(tree_build)
    
    ###--------multi-process------------------
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


