#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2017, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    1. Identify SSRs of given sequences
    2. Extract flanking sequences of identified SSRs
    3. Design primers for each SSR regions
    4. Check amplicon status of each primer

Program requirement:
    1. MISA
    2. primer3
    3. eprimer32
    4. primersearch

Fatsa seq:

>seq_id
ACGCTACGACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTGCAGAGAGTGAGATG
ACGCATTGAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTACGAAATAGGAA

SSR file (output by MISA)

ID      SSR nr. SSR type        SSR     size    start   end
comp1_c0_seq1   1       p1      (A)15   15      495     509
comp16_c0_seq1  1       p2      (GA)11  22      1       22
comp16_c0_seq2  1       p2      (GA)11  22      1       22
comp18_c0_seq1  1       p1      (T)11   11      8       18
comp24_c0_seq1  1       p5      (TAGCC)5        25      558     582

'''

import sys
import os
import math
#import subprocess as sub
#from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

def outputMISA(misa="misa.pl"):
    #if os.path.exists(misa):
    #    return

    misa_fh = open(misa,'w')
    print >>misa_fh, """#!/usr/bin/perl -w
# Author: Thomas Thiel
# Program name: misa.pl

###_______________________________________________________________________________
###
###Program name: misa.pl
###Author:       Thomas Thiel
###Release date: 14/12/01 (version 1.0)
###
###_______________________________________________________________________________
###
## _______________________________________________________________________________
##
## DESCRIPTION: Tool for the identification and localization of
##              (I)  perfect microsatellites as well as
##              (II) compound microsatellites (two individual microsatellites,
##                   disrupted by a certain number of bases)
##
## SYNTAX:   misa.pl <FASTA file>
##
##    <FASTAfile>    Single file in FASTA format containing the sequence(s).
##
##    In order to specify the search criteria, an additional file containing
##    the microsatellite search parameters is required named "misa.ini", which
##    has the following structure:
##      (a) Following a text string beginning with 'def', pairs of numbers are
##          expected, whereas the first number defines the unit size and the
##          second number the lower threshold of repeats for that specific unit.
##      (b) Following a text string beginning with 'int' a single number defines
##          the maximal number of bases between two adjacent microsatellites in
##          order to specify the compound microsatellite type.
##    Example:
##      definition(unit_size,min_repeats):          1-10 2-6 3-5 4-5 5-5 6-5
##      interruptions(max_difference_for_2_SSRs):   100
##
## EXAMPLE: misa.pl seqs.fasta
##
## _______________________________________________________________________________
##



# Check for arguments. If none display syntax #

if (@ARGV == 0)
  {
  open (IN,"<$0");
  while (<IN>) {if (/^\#\# (.*)/) {$message .= "$1\\n"}};
  close (IN);
  die $message;
  };

# Check if help is required #

if ($ARGV[0] =~ /-help/i)
  {
  open (IN,"<$0");
  while (<IN>) {if (/^\#\#\#(.*)/) {$message .= "$1\\n"}};
  close (IN);
  die $message;
  };

# Open FASTA file #

open (IN,"<$ARGV[0]") || die ("\\nError: FASTA file doesn't exist !\\n\\n");
open (OUT,">$ARGV[0].misa");
print OUT "ID\\tSSR nr.\\tSSR type\\tSSR\\tsize\\tstart\\tend\\n";

# Reading arguments #

open (SPECS,"misa.ini") || die ("\\nError: Specifications file doesn't exist !\\n\\n");
my %typrep;
my $amb = 0;
while (<SPECS>)
   {
   %typrep = $1 =~ /(\d+)/gi if (/^def\S*\s+(.*)/i);
   if (/^int\S*\s+(\d+)/i) {$amb = $1}
   };
my @typ = sort { $a <=> $b } keys %typrep;


# CORE

$/ = ">";
my $max_repeats = 1; #count repeats
my $min_repeats = 1000; #count repeats
my (%count_motif,%count_class); #count
my ($number_sequences,$size_sequences,%ssr_containing_seqs); #stores number and size of all sequences examined
my $ssr_in_compound = 0;
my ($id,$seq);
while (<IN>)
  {
  next unless (($id,$seq) = /(.*?)\\n(.*)/s);
  my ($nr,%start,@order,%end,%motif,%repeats); # store info of all SSRs from each sequence
  $seq =~ s/[\d\s>]//g; #remove digits, spaces, line breaks,...
  $id =~ s/^\s*//g; $id =~ s/\s*$//g;$id =~ s/\s/_/g; #replace whitespace with "_"
  $number_sequences++;
  $size_sequences += length $seq;
  for ($i=0; $i < scalar(@typ); $i++) #check each motif class
    {
    my $motiflen = $typ[$i];
    my $minreps = $typrep{$typ[$i]} - 1;
    if ($min_repeats > $typrep{$typ[$i]}) {$min_repeats = $typrep{$typ[$i]}}; #count repeats
    my $search = "(([acgt]{$motiflen})\\\\2{$minreps,})";
    while ( $seq =~ /$search/ig ) #scan whole sequence for that class
      {
      my $motif = uc $2;
      my $redundant; #reject false type motifs [e.g. (TT)6 or (ACAC)5]
      for ($j = $motiflen - 1; $j > 0; $j--)
        {
        my $redmotif = "([ACGT]{$j})\\\\1{".($motiflen/$j-1)."}";
        $redundant = 1 if ( $motif =~ /$redmotif/ )
        };
      next if $redundant;
      $motif{++$nr} = $motif;
      my $ssr = uc $1;
      $repeats{$nr} = length($ssr) / $motiflen;
      $end{$nr} = pos($seq);
      $start{$nr} = $end{$nr} - length($ssr) + 1;
      # count repeats
      $count_motifs{$motif{$nr}}++; #counts occurrence of individual motifs
      $motif{$nr}->{$repeats{$nr}}++; #counts occurrence of specific SSR in its appearing repeat
      $count_class{$typ[$i]}++; #counts occurrence in each motif class
      if ($max_repeats < $repeats{$nr}) {$max_repeats = $repeats{$nr}};
      };
    };
  next if (!$nr); #no SSRs
  $ssr_containing_seqs{$nr}++;
  @order = sort { $start{$a} <=> $start{$b} } keys %start; #put SSRs in right order
  $i = 0;
  my $count_seq; #counts
  my ($start,$end,$ssrseq,$ssrtype,$size);
  while ($i < $nr)
    {
    my $space = $amb + 1;
    if (!$order[$i+1]) #last or only SSR
      {
      $count_seq++;
      my $motiflen = length ($motif{$order[$i]});
      $ssrtype = "p".$motiflen;
      $ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}";
      $start = $start{$order[$i]}; $end = $end{$order[$i++]};
      next
      };
    if (($start{$order[$i+1]} - $end{$order[$i]}) > $space)
      {
      $count_seq++;
      my $motiflen = length ($motif{$order[$i]});
      $ssrtype = "p".$motiflen;
      $ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}";
      $start = $start{$order[$i]}; $end = $end{$order[$i++]};
      next
      };
    my ($interssr);
    if (($start{$order[$i+1]} - $end{$order[$i]}) < 1)
      {
      $count_seq++; $ssr_in_compound++;
      $ssrtype = 'c*';
      $ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}($motif{$order[$i+1]})$repeats{$order[$i+1]}*";
      $start = $start{$order[$i]}; $end = $end{$order[$i+1]}
      }
    else
      {
      $count_seq++; $ssr_in_compound++;
      $interssr = lc substr($seq,$end{$order[$i]},($start{$order[$i+1]} - $end{$order[$i]}) - 1);
      $ssrtype = 'c';
      $ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}$interssr($motif{$order[$i+1]})$repeats{$order[$i+1]}";
      $start = $start{$order[$i]};  $end = $end{$order[$i+1]};
      #$space -= length $interssr
      };
    while ($order[++$i + 1] and (($start{$order[$i+1]} - $end{$order[$i]}) <= $space))
      {
      if (($start{$order[$i+1]} - $end{$order[$i]}) < 1)
        {
        $ssr_in_compound++;
        $ssrseq .= "($motif{$order[$i+1]})$repeats{$order[$i+1]}*";
        $ssrtype = 'c*';
        $end = $end{$order[$i+1]}
        }
      else
        {
        $ssr_in_compound++;
        $interssr = lc substr($seq,$end{$order[$i]},($start{$order[$i+1]} - $end{$order[$i]}) - 1);
        $ssrseq .= "$interssr($motif{$order[$i+1]})$repeats{$order[$i+1]}";
        $end = $end{$order[$i+1]};
        #$space -= length $interssr
        }
      };
    $i++;
    }
  continue
    {
    print OUT "$id\\t$count_seq\\t$ssrtype\\t$ssrseq\\t",($end - $start + 1),"\\t$start\\t$end\\n"
    };
  };

close (OUT);
open (OUT,">$ARGV[0].statistics");

# INFO

# Specifications

print OUT "Specifications\\n==============\\n\\nSequence source file: \\"$ARGV[0]\\"\\n\\nDefinement of microsatellites (unit size / minimum number of repeats):\\n";
for ($i = 0; $i < scalar (@typ); $i++) {print OUT "($typ[$i]/$typrep{$typ[$i]}) "};print OUT "\\n";
if ($amb > 0) {print OUT "\\nMaximal number of bases interrupting 2 SSRs in a compound microsatellite:  $amb\\n"};
print OUT "\\n\\n\\n";

# OCCURRENCE OF SSRs

#small calculations
my @ssr_containing_seqs = values %ssr_containing_seqs;
my $ssr_containing_seqs = 0;
for ($i = 0; $i < scalar (@ssr_containing_seqs); $i++) {$ssr_containing_seqs += $ssr_containing_seqs[$i]};
my @count_motifs = sort {length ($a) <=> length ($b) || $a cmp $b } keys %count_motifs;
my @count_class = sort { $a <=> $b } keys %count_class;
for ($i = 0; $i < scalar (@count_class); $i++) {$total += $count_class{$count_class[$i]}};


# Overview

print OUT "RESULTS OF MICROSATELLITE SEARCH\\n================================\\n\\n";
print OUT "Total number of sequences examined:              $number_sequences\\n";
print OUT "Total size of examined sequences (bp):           $size_sequences\\n";
print OUT "Total number of identified SSRs:                 $total\\n";
print OUT "Number of SSR containing sequences:              $ssr_containing_seqs\\n";
print OUT "Number of sequences containing more than 1 SSR:  ",$ssr_containing_seqs - ($ssr_containing_seqs{1} || 0),"\\n";
print OUT "Number of SSRs present in compound formation:    $ssr_in_compound\\n\\n\\n";


# Frequency of SSR classes

print OUT "Distribution to different repeat type classes\\n---------------------------------------------\\n\\n";
print OUT "Unit size\\tNumber of SSRs\\n";
my $total = undef;
for ($i = 0; $i < scalar (@count_class); $i++) {print OUT "$count_class[$i]\\t$count_class{$count_class[$i]}\\n"};
print OUT "\\n";


# Frequency of SSRs: per motif and number of repeats

print OUT "Frequency of identified SSR motifs\\n----------------------------------\\n\\nRepeats";
for ($i = $min_repeats;$i <= $max_repeats; $i++) {print OUT "\\t$i"};
print OUT "\\ttotal\\n";
for ($i = 0; $i < scalar (@count_motifs); $i++)
  {
  my $typ = length ($count_motifs[$i]);
  print OUT $count_motifs[$i];
  for ($j = $min_repeats; $j <= $max_repeats; $j++)
    {
    if ($j < $typrep{$typ}) {print OUT "\\t-";next};
    if ($count_motifs[$i]->{$j}) {print OUT "\\t$count_motifs[$i]->{$j}"} else {print OUT "\\t"};
    };
  print OUT "\\t$count_motifs{$count_motifs[$i]}\\n";
  };
print OUT "\\n";


# Frequency of SSRs: summarizing redundant and reverse motifs

# Eliminates %count_motifs !
print OUT "Frequency of classified repeat types (considering sequence complementary)\\n-------------------------------------------------------------------------\\n\\nRepeats";
my (%red_rev,@red_rev); # groups
for ($i = 0; $i < scalar (@count_motifs); $i++)
  {
  next if ($count_motifs{$count_motifs[$i]} eq 'X');
  my (%group,@group,$red_rev); # store redundant/reverse motifs
  my $reverse_motif = $actual_motif = $actual_motif_a = $count_motifs[$i];
  $reverse_motif =~ tr/ACGT/TGCA/;
  $reverse_motif = reverse $reverse_motif;
  my $reverse_motif_a = $reverse_motif;
  for ($j = 0; $j < length ($count_motifs[$i]); $j++)
    {
    if ($count_motifs{$actual_motif}) {$group{$actual_motif} = "1"; $count_motifs{$actual_motif}='X'};
    if ($count_motifs{$reverse_motif}) {$group{$reverse_motif} = "1"; $count_motifs{$reverse_motif}='X'};
    $actual_motif =~ s/(.)(.*)/$2$1/;
    $reverse_motif =~ s/(.)(.*)/$2$1/;
    $actual_motif_a = $actual_motif if ($actual_motif lt $actual_motif_a);
    $reverse_motif_a = $reverse_motif if ($reverse_motif lt $reverse_motif_a)
    };
  if ($actual_motif_a lt $reverse_motif_a) {$red_rev = "$actual_motif_a/$reverse_motif_a"}
  else {$red_rev = "$reverse_motif_a/$actual_motif_a"}; # group name
  $red_rev{$red_rev}++;
  @group = keys %group;
  for ($j = 0; $j < scalar (@group); $j++)
    {
    for ($k = $min_repeats; $k <= $max_repeats; $k++)
      {
      if ($group[$j]->{$k}) {$red_rev->{"total"} += $group[$j]->{$k};$red_rev->{$k} += $group[$j]->{$k}}
      }
    }
  };
for ($i = $min_repeats; $i <= $max_repeats; $i++) {print OUT "\\t$i"};
print OUT "\\ttotal\\n";
@red_rev = sort {length ($a) <=> length ($b) || $a cmp $b } keys %red_rev;
for ($i = 0; $i < scalar (@red_rev); $i++)
  {
  my $typ = (length ($red_rev[$i])-1)/2;
  print OUT $red_rev[$i];
  for ($j = $min_repeats; $j <= $max_repeats; $j++)
    {
    if ($j < $typrep{$typ}) {print OUT "\\t-";next};
    if ($red_rev[$i]->{$j}) {print OUT "\\t",$red_rev[$i]->{$j}}
    else {print OUT "\\t"}
    };
  print OUT "\\t",$red_rev[$i]->{"total"},"\\n";
  };

"""
    misa_fh.close()
#----outputMISA---------------------------------------------------------------------

def output_misa_ini(misa_ini="misa.ini"):
    #if os.path.exists(misa_ini):
    #    return

    misa_ini_fh = open(misa_ini,'w')
    print >>misa_ini_fh, """definition(unit_size,min_repeats):                   1-10 2-6 3-5 4-5 5-5 6-5
interruptions(max_difference_between_2_SSRs):        0"""
    misa_ini_fh.close()

#--------------END of output_misa_ini------------------

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
#-------------------------------------------------------

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
        metavar="FILEIN", help="FASTA file used for SSR identification. Blanks in FASTA id will be transferred into '_'.")
    parser.add_option("-d", "--databse-file", dest="database",
        help="FASTA file containing multiple speceis info")
    parser.add_option("-s", "--ssr-file", dest="ssr",
        help="Optional. SSR file generated by misa.pl. If not given, SSR will be first identified.")
    parser.add_option("-S", "--amplicon-size", dest="amplicon_size",
        default = "200,100,280", 
        help="Amplicon size. Default <200,100,280> represents <optampliconsize,minampliconsize,maxampliconsize>. It will affect both SSR region extracttion and primer design.")
    parser.add_option("-f", "--flank", dest="flank",
        default = "13,250", 
        help="Length of flank regions to be extracted along each SSR. Default <13,250> meaning at least 13 nt and at most 250 nt will be extracted.")
    parser.add_option("-P", "--primer-size", dest="primer_size",
        default = "20,15,25", 
        help="Default <20,15,25> represents <optsize,minsize,maxsize> respectively")
    parser.add_option("-T", "--tm", dest="primer_tm",
        default = "50,45,55", 
        help="Default <50,45,55> represents <opttm,mintm,maxtm> respectively")
    parser.add_option("-p", "--output-prefix", dest="op",
        help="Output prefix. Optional")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readmisaSSR(misa_ssr):
    """
    ssrDict = {'scaffold': 
                    [
                        {
                            ssr_type: "p1",
                            ssr: "(A)15",
                            size: 15,
                            start: 495,
                            end: 509
                        },
                        {
                            ssr_type: "p1",
                            ssr: "(A)15",
                            size: 15,
                            start: 495,
                            end: 509
                        },
                    ]
                }
    """
    header = 1
    ssrDict = {}
    for line in open(misa_ssr):
        if header:
            header -= 1
            continue
        #-------------------
        lineL = line.split('\t')
        seqId = lineL[0]
        ssr_type = lineL[2]
        ssr      = lineL[3]
        size     = int(lineL[4])
        # Original: 1 -start, both included
        # Transfer to : 0-start, start but not end included
        start    = int(lineL[5]) - 1
        end      = int(lineL[6]) 
        if seqId not in ssrDict:
            ssrDict[seqId] = []
        subD = {}
        subD["ssr_type"] = ssr_type
        subD["ssr"] = ssr
        subD["size"] = size
        subD["start"] = start
        subD["end"] = end
        ssrDict[seqId].append(subD)
    #--------------------------------
    return ssrDict
#-------------END readmisaSSR--------------

#"""
# amplicon size between 100
# and 280 bp; minimum, optimum, and maximum annealing
# temperature (TM) of 45, 50, and 55, respectively, minimum,
# optimum, and maximum primer size of 15, 20, and 25 bp,
# respectively.
#"""

def extractSeq(output_seqfh,key,seq,ssrL,flankMin,flankMax, am_size_min, am_size_max):
    """
    ssrL = 
        [
            {
                ssr_type: "p1",
                ssr: "(A)15",
                size: 15,
                start: 495,
                end: 509
            },
            {
                ssr_type: "p1",
                ssr: "(A)15",
                size: 15,
                start: 495,
                end: 509
            },
        ]
    """
    len_seq = len(seq)
    for ssrDsub in ssrL:
        start = ssrDsub['start']
        end = ssrDsub['end']
        ssr_type = ssrDsub['ssr_type']
        ssr = ssrDsub['ssr']
        size = ssrDsub['size']
        
        # Deplete very short sequences
        if len_seq < am_size_min:
            continue
        # Delete too small upstream and down stream sequences
        if start < flankMin or end + flankMin > len_seq:
            continue
        # Re-asign length
        up_potential = start
        end_potential = len_seq - end
        half_am_size_max = am_size_max/2
        if up_potential >= half_am_size_max and end_potential >=half_am_size_max:
            flank_max = half_am_size_max
        elif up_potential <= half_am_size_max and end_potential >= half_am_size_max:
            flank_max = am_size_max - up_potential
        elif up_potential >= half_am_size_max and end_potential <= half_am_size_max:
            flank_max = am_size_max - end_potential
        #--------------------------------------

        ssr_seq = seq[start:end].lower()
        up_start = start - flankMax
        if up_start < 0:
            up_start = 0
        up_end = start
        #print >>sys.stderr, up_start,up_end
        ssr_up_seq = seq[up_start:up_end]
        dw_start = end
        dw_end   = dw_start + flankMax
        if dw_end > len_seq:
            dw_end = len_seq
        ssr_dw_seq = seq[dw_start:dw_end]
        #-----------------------------------
        print >>output_seqfh, ">{} {} {}".format(key, up_end-up_start-1,dw_start-up_start+1)
        #print seq
        #print "#"
        #print ssr_up_seq
        #print "##"
        #print ssr_seq
        #print "###"
        #print ssr_dw_seq
        print >>output_seqfh, ''.join([ssr_up_seq, ssr_seq, ssr_dw_seq])
#-----------extractSeq--------------------

def extractFlankSeq(file, output_seqFile, ssrDict, flankMin, flankMax, am_size_min, am_size_max):
    """
    ssrDict = {'scaffold': 
                    [
                        {
                            ssr_type: "p1",
                            ssr: "(A)15",
                            size: 15,
                            start: 495,
                            end: 509
                        },
                        {
                            ssr_type: "p1",
                            ssr: "(A)15",
                            size: 15,
                            start: 495,
                            end: 509
                        },
                    ]
                }
    """
    output_seqfh = open(output_seqFile, 'w')
    seqL = []
    for line in open(file):
        line = line.strip()
        if line[0]== '>':
            if seqL:
                seq = ''.join(seqL)
                ssrL = ssrDict.get(key,[])
                if not ssrL:
                    key2 = key.replace(' ','_')
                    ssrL = ssrDict.get(key2,[])
                if ssrL:
                    extractSeq(output_seqfh,key,seq,ssrL,flankMin,flankMax, am_size_min, am_size_max)
            key = line[1:]
            seqL = []
        else:
            seqL.append(line)
        #--------------------------------------
    #-------------------------------------
    if seqL:
        seq  = ''.join(seqL)
        ssrL = ssrDict.get(key,[])
        if not ssrL:
            key2 = key.replace(' ','_')
            ssrL = ssrDict.get(key2,[])
        if ssrL:
            extractSeq(output_seqfh, key,seq,ssrL,flankMin,flankMax, am_size_min, am_size_max)
    #---------------------------------------
    output_seqfh.close()
#-------------------------------------------------------


def transferPrimer3OutputToPrimersearchInput(primer3, primersearch_fh):
    fh = open(primer3)
    #fh_out = open(primersearch,'w')
    line = fh.readline()
    line = fh.readline()
    assert line.find("EPRIMER32 RESULTS FOR") != -1, "Wrong format" + primer3
    seq_name = line.strip()[24:]
    
    #primerL = []
    count = 1
    for line in fh:
        if line.find("FORWARD PRIMER") != -1:
            forward = line.strip().split()[-1]
        if line.find("REVERSE PRIMER") != -1:
            reverse = line.strip().split()[-1]
            #tmpL = [seq_name+'@'+str(count), forward, reverse]
            print >>primersearch_fh, "{}@{}\t{}\t{}".format(seq_name,count,forward, reverse)
            #primerL.add(tmpL)
            count += 1
    #------------------------------------------
    fh.close()
    #fh_out.close()
#---------------------------------------------------

def deisgnPrimerForEachSeq(output_seq, primersearch_fh, optsize, minsize, maxsize, opttm, mintm, maxtm, am_size_opt, am_size_min, am_size_max):
    """
    Single line FASTA file
    """
    count = 1
    for line in open(output_seq):
        if line[0] == '>':
            id,start,end = line.strip().split(' ')
            fileid = str(count)
            count += 1
        else:
            output_file = output_seq + fileid + '.fa'
            output_file_fh = open(output_file, 'w')
            print >>output_file_fh, "{}\n{}".format(id, line.strip())
            output_file_fh.close()
            cmd = ["eprimer32 -sequence", output_file, '-outfile', output_file+'.primer', 
                    "-targetregion", start+','+end, '-optsize', optsize, "-numreturn 3",
                    '-minsize', minsize, '-maxsize', maxsize, "-opttm", opttm, "-mintm", mintm, "-maxtm", maxtm,
                    "-psizeopt", str(am_size_opt), "-prange", str(am_size_min)+'-'+str(am_size_max),
                    '-die -auto']
            cmd = ' '.join(cmd)
            #print >>sys.stderr, cmd
            #p = sub.Popen(cmd, stdout=sub.PIPE,stderr=sub.PIPE)
            #output, errors = p.communicate()
            #print >>sys.stderr,output
            #print >>sys.stderr,errors
            #os.system("pwd")
            if os.system(cmd):
                print >>sys.stderr, cmd + ' Wrong'
                sys.exit(1)
            #--------------------------------------------------
            transferPrimer3OutputToPrimersearchInput(output_file+'.primer', primersearch_fh)
    #-------------END for-----------------------------
#-----------------------------------------------------


def readInPrimerSearch(file, primerDict, species=''):
    '''
    Primer name SRR037890___comp24_c0_seq1@1
    Amplimer 1
        Sequence: SRR037890___comp24_c0_seq1  
        
        TATTTTCCTATGTTGCTACC hits forward strand at 433 with 0 mismatches
        ACTATTAGCTGTAAAGCAAA hits reverse strand at [23] with 0 mismatches
        Amplimer length: 200 bp


    primerDict = {
        primerName : {

            # Record amplicon info
            # amplicon_start, amplicon_end: 0-started, second excluded
            'species': [
                [targetseq1, AmpliconSize, (amplicon_start, amplicon_end), (forward_mismatch, reverse_mismatch)]
                [targetseq2, AmpliconSize, (amplicon_start, amplicon_end), (forward_mismatch, reverse_mismatch)]
            ],

            # Summary amplicon info
            # Compute number of amplicons for each size of each species
            'sta':{
                species1: {
                    200: 1,
                    300, 1
                },
                species2: {
                    200: 2,
                    300, 1
                }
            }
        }
    }
    '''
    
    seqD = {}
    for line in open(file):
        if line.find("Primer name") == 0:
            primer_name = line[12:-1]
            if primer_name not in primerDict:
                primerDict[primer_name] = {}
                primerDict[primer_name]['sta'] = {}
        elif line.find("Sequence:") != -1:
            targetSeq = line.strip().split()[1]
            seqD[targetSeq] = 1
            if not species:
                species = targetSeq.split('___')[0]
            seq_name = targetSeq
            if species not in primerDict[primer_name]:
                primerDict[primer_name][species] =[]
                primerDict[primer_name]['sta'][species] ={}
        elif line.find("hits forward strand at")  != -1:
            lineL = line.strip().split()
            forward_mismatch = int(lineL[-2])
            amplicon_start = int(lineL[-4])-1
        elif line.find("hits reverse strand at") != -1:
            lineL = line.strip().split()
            reverse_mismatch = int(lineL[-2])
            amplicon_end = (-1) * int(lineL[-4][1:-1])+1
        elif line.find("Amplimer length") != -1:
            ampliconSize = int(line.strip().split()[2])
            primerDict[primer_name][species].append([seq_name, ampliconSize, (amplicon_start, amplicon_end), (forward_mismatch, reverse_mismatch)])
            primerDict[primer_name]['sta'][species][ampliconSize] = primerDict[primer_name]['sta'][species].get(ampliconSize,0)+1
        #----------------------------------------------------------          
        #print >>sys.stderr,primerDict
    #------------------END for---------------------
    return seqD
#-------------------------------------------------------

def gamma(x, t, k=2):
    b = 1.0 / t
    x = x * k
    c = 0.0
    
    for i in range(0, k):
        c += (math.exp(-b * x) * (b * x) ** i) / math.factorial(i)

    return c
#-------END gamma-------------------
def readInGelFile(file):
    '''
    Molecular weight or amplicon size

    SIZE    samp1   samp2   samp3   samp4
    size1   0   0   0   0
    size2   200   0   0   0
    size3   200   200   0   200
    size4   0   0   200   0
    '''
    
    header = 1
    aDict = {}
    for line in open(file):
        if header:
            sampleL = line.split()
            for sample in sampleL[1:]:
                aDict[sample] = []
            header -= 1
            continue
        #_-----------------------
        lineL = line.split()
        key = float(lineL[0])
        for sample, concentration in zip(sampleL[1:], lineL[1:]):
            concentration = int(concentration)
            if concentration:
                aDict[sample].append([key,concentration])
        #----------------------------------------------
    return aDict,sampleL[1:]
#----------------------------


def plot(band_matrix, xpositionS, xvalueS, markerPositionS, markerSizeS):

    output = band_matrix + '.r'
    output_fh = open(output, 'w')
    print >>output_fh, """
usePackage <- function(p) {{
	    if (!is.element(p, installed.packages()[,1]))
			        install.packages(p, dep = TRUE)
 	   require(p, character.only = TRUE)
}}

usePackage("ggplot2")
usePackage("reshape2")


data <- read.table(file="{file}", sep="\t", header=T, row.names=1,
	check.names=F, quote="")


data$id <- rownames(data)
idlevel <- as.vector(rownames(data))

idlevel <- rev(idlevel)

data.m <- melt(data, c("id"))

data.m$id <- factor(data.m$id, levels=idlevel, ordered=T)

p <- ggplot(data=data.m, aes(x=variable, y=id)) + 	geom_tile(aes(fill=value)) 

midpoint = 0

p <- p + scale_fill_gradient2(low="black", mid="grey",
	high="white", midpoint=midpoint, name=" ",
	na.value="grey")

p <- p + theme(axis.ticks=element_blank()) + theme_bw() + 
	theme(panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), panel.border = element_blank()) + xlab("") + 
	ylab("") + labs(title="")

p <- p + theme(axis.text.x=element_text(angle=45,hjust=0, vjust=0))

none='none'
legend_pos_par <- none

p <- p + theme(legend.position=legend_pos_par)


xtics_pos <- {xpositionS}
xtics_value <- {xvalueS}

p <- p + scale_x_discrete(breaks=xtics_pos, labels=xtics_value, position="top")

ytics_pos <- {markerPositionS}
ytics_value <- {markerSizeS}

p <- p + scale_y_discrete(breaks=ytics_pos, labels=ytics_value, position="right")

p <- p + theme(text=element_text(size=14))


ggsave(p, filename="{file}.png", dpi=300, width=10, height=20, units=c("cm"))
""".format(file=band_matrix, xvalueS=xvalueS, xpositionS=xpositionS, markerSizeS=markerSizeS, 
        markerPositionS=markerPositionS)

    output_fh.close()
    if os.system("Rscript "+output):
        print >>sys.stderr, "Wrong in running plot"

#----------------------plot----------------------------------------------

def gel_main(file, gel_concentration=2, voltage=20, time=40):
    '''
    gel_concentration: Default <2> represents 2%
    voltage: 20 v
    time: 40 minutes
    '''
    optimum_DNA_length = 2000 / gel_concentration ** 3

    sampleD, nameL = readInGelFile(file)

    # Add ladder
    nameL.append("Markers")
    sampleD["Markers"] = [(50,200), (100,200),(150,200),(200,200),(250,200),(300,200),(350,200),(400,200),(500,200),(1000,200)]

    len_sampleD = len(sampleD)


    # Number of lanes
    lane_count = len_sampleD
    # Width and height of dingle lane
    lane_width = 30
    lane_height = 3
    # Intervals between neighboring lanes
    lane_interval = 6
    gel_border  = 12

    gel_width = 2 * gel_border + lane_count * lane_width + (lane_count-1) * lane_interval

    #gel_height = gel_width

    # Generate empty lanes
    # Every 1 unit from up_border to all height
    strandD = {}


    # Loading samples
    # position: represents start position, loaing hole
    for name in nameL:
        strandD[name] = []
        dnaL = sampleD[name]
        for dna in dnaL:
            band = {"size": float(dna[0]), 'conc': dna[1]*1.0/lane_height, 'position': gel_border+5}
            strandD[name].append(band)
    
    # Run 
    # Move loaded DNA down the gel at a rate dependent on voltage, DNA length and agarose concentration.
    time = time
    voltage = voltage
    max_dist = 0.25 * time * voltage

    maxposition = 0
    for name in nameL:
        for bandD in strandD[name]:
            g = gamma(bandD['size']/20, int(optimum_DNA_length/20))   
            bandD['position'] += max_dist * g
            if maxposition < bandD['position']:
                maxposition = bandD['position']
            #print >>sys.stderr, bandD
    #-----------------------------------------------------
    maxposition = int(maxposition+gel_border+1)
    # Determines where in the concentration of DNA in every part of the gel
    # Generate an empty lane, O represents no band

    markerSize     = [str(int(bandD["size"])) for bandD in strandD["Markers"]]
    markerSize.reverse()
    markerSizeS = "c("+','.join(["'"+i+"'" for i in markerSize])+ ")"
    markerPosition = ['CT'+str(int(bandD["position"]+0.5)) for bandD in strandD["Markers"]]
    markerPosition.reverse()
    markerPositionS = "c("+','.join(["'"+i+"'" for i in markerPosition])+ ")"
    #print >>sys.stderr, markerSizeS
    #print >>sys.stderr, markerPositionS
    #laneL = [[-1 for j in range(maxposition)] for i in range(len_sampleD)]

    laneL = []
    for i in range(len_sampleD):
        tmpL = [-1000 for j in range(maxposition)]
        tmpL[gel_border] = -10
        tmpL[gel_border-1] = 0
        tmpL[gel_border-2] = -10
        laneL.append(tmpL)
    #-------------------------------------
    band_count = maxposition

    for i in range(len_sampleD):
        name = nameL[i]
        for bandD in strandD[name]:
            for y in range(lane_height-2):
                pos = int(bandD['position'])+y
                if pos < band_count - 4:
                    laneL[i][pos-1] += 0.12 * bandD['conc'] * bandD['size']
                    laneL[i][pos]   += 0.2  * bandD['conc'] * bandD['size']
                    laneL[i][pos+1] += 0.12 * bandD['conc'] * bandD['size']
                #-----Blur edges-----------------------------
            #---------------------------
        #_--------------------------
    #-------------------------------
    max_value = max([max(i) for i in laneL])
    min_value = -1 * max_value
    for i in range(len_sampleD):
        for j in range(maxposition):
            if laneL[i][j] == -1000:
                laneL[i][j] = min_value
    #print laneL    
    # Expose
    # print "ID\t{}".format("\t".join(nameL))
    newCol = len(nameL) * 8 + 2
    xposition = ['ct'+str(i) for i in range(3, newCol,8)]

    xpositionS = "c("+','.join(["'"+i+"'" for i in xposition])+ ")"
    xvalueS = "c("+','.join(["'"+i+"'" for i in nameL])+ ")"
    #print >>sys.stderr, xpositionS
    #print >>sys.stderr, '\t'.join(nameL)

    # Plot data

    band_matrix = file + '.virtualGel.data'
    band_matrix_fh = open(band_matrix, 'w')
    print >>band_matrix_fh, "ID\t{}".format('\t'.join(['ct'+str(i) for i in range(newCol)]))
    for i in range(maxposition):
        tmpL = ["CT"+str(i)]
        tmpL.append(str(min_value))
        tmpL.append(str(min_value))
        for lane in laneL:
            # Blur edges, 6 main lane, two margin lane
            if lane[i] > 0:
                tmpL.append(str(0.8*lane[i]))
                tmpL.append(str(0.9*lane[i]))
            else:
                tmpL.append(str(lane[i]))
                tmpL.append(str(lane[i]))
            tmpL.append(str(lane[i]))
            tmpL.append(str(lane[i]))

            if lane[i] > 0:
                tmpL.append(str(0.9*lane[i]))
                tmpL.append(str(0.8*lane[i]))
            else:
                tmpL.append(str(lane[i]))
                tmpL.append(str(lane[i]))

            tmpL.append(str(min_value))
            tmpL.append(str(min_value))
        print >>band_matrix_fh, "\t".join(tmpL)
    band_matrix_fh.close()
    # Generate R script for plot 
    
    plot(band_matrix, xpositionS, xvalueS, markerPositionS, markerSizeS)
#---------END main--------------------------------

def plotamplicon(finalOutput, primer_id, ampliconsizeL, targetSpeL):
    '''
    ampliconsizeL = [
        Query_spe:100,200,300,
        target_spe1:100,
        target_spe2:100,
        target_spe3:100,
    ]
    '''
    gel_plot_data = finalOutput+primer_id+'.gel'
    gel_plot_data_fh = open(gel_plot_data, 'w')
    
    targetSpeL.insert(0, 'Query_spe')
    
    #sizeL = [i.split(':') for i in ampliconsizeL]
    sizeD = {}
    for each in ampliconsizeL:
        spe, size = each.rsplit(':',1)
        sizeL = [int(i) for i in size.split(',')]
        for i in sizeL:
            if i not in sizeD:
                sizeD[i] = {}
            sizeD[i][spe] = '200'
    print >>gel_plot_data_fh, "ID\t{}".format('\t'.join(targetSpeL))
    
    sizeL = sizeD.keys()
    sizeL.sort()

    for size in sizeL:
        tmpL = [sizeD[size].get(spe,'0') for spe in targetSpeL]
        print >>gel_plot_data_fh, "{}\t{}".format(size, '\t'.join(tmpL))
    gel_plot_data_fh.close()
    
    gel_main(gel_plot_data)

    return gel_plot_data+'.virtualGel.data.png'
#----------------------------------------------

def outputampliconSeq(finalOutput, primer_id, speD, query_fastaD, matchSpeL, target_fastaD):
    ampliconSeqFile = finalOutput + '.' + primer_id + '.fa'
    ampliconSeqFile_fh = open(ampliconSeqFile, 'w')

    querySeqL = speD['Query_spe']
    #queryampliconL = []
    count = 1
    for eachQueryMatch in querySeqL:
        id = eachQueryMatch[0]
        start,end = eachQueryMatch[2]
        print >>ampliconSeqFile_fh, \
            ">%s %d\n%s" % ('Query_spe', count, query_fastaD[id][start:end])
        count += 1
        #queryampliconL.append(query_fastaD[id][start:end])
    #queryAmplicon = 'Query_spe'+':'+','.join(queryampliconL)
    #-----------------------------------------------------------
    #ampliconL = [queryAmplicon]
    for otherspe in matchSpeL:
        count = 1
        targetSeqL = speD[otherspe]
        #tmpL = []
        for eachQueryMatch in targetSeqL:
            id = eachQueryMatch[0]
            start,end = eachQueryMatch[2]
            #tmpL.append(target_fastaD[id][start:end])
            print >>ampliconSeqFile_fh, \
                ">%s %d\n%s" % (otherspe, count, target_fastaD[id][start:end])
            count += 1
        #targetAmplicon = otherspe+':'+','.join(tmpL)
        #ampliconL.append(targetAmplicon)
    ampliconSeqFile_fh.close()
    return ampliconSeqFile
#---------------------------------------------------------------------

def tabulizePrimerSearch(primerDict, primerSeqD, query_fastaD, target_fastaD, targetSpeL, finalOutput):
    '''
    '''
    targetSpeNum = len(targetSpeL)
    finalOutput_fh = open(finalOutput, 'w')
    primeridL = primerDict.keys()
    primeridL.sort(key = lambda x: len(primerDict[x]['sta']), reverse=True)
    headerL = ["Primer ID", "Match species", "Match No", "differentAmpliconCnt", "Amplicon size", "Forward primer", "Reverse Primer", "Amplicon sequence"]
    for primer_id in primeridL:
        speD = primerDict[primer_id]
        staD = speD['sta']
        #clone_spe = len(staD)
        #target = speD['Query_spe']
        matchSpeL = speD.keys()
        matchSpeL.remove("sta")
        matchSpeL.remove('Query_spe')
        clone_spe = len(matchSpeL)
        matchSpeL.sort()
        
        ampliconSizeL = [subspe+':'+','.join([str(i) for i in staD[subspe].keys()]) for subspe in matchSpeL]
        queryAmpliconSize = 'Query_spe'+':'+ ','.join([str(i) for i in staD['Query_spe'].keys()])
        ampliconSizeL.insert(0, queryAmpliconSize)
        
        plot_name = plotamplicon(finalOutput, primer_id, ampliconSizeL[:], targetSpeL[:])

        differentAmpliconCnt = len(set(ampliconSizeL))
        
        ampliconsize =';'.join(ampliconSizeL)+'|'+plot_name
        #print >>sys.stderr, primerSeqD[primer_id]
        forward_primer, reverse_primer = primerSeqD[primer_id].split()
        
        #----------------------------------------------------------------------
        ampliconFile = outputampliconSeq(finalOutput, primer_id, speD, query_fastaD, matchSpeL, target_fastaD)
        #querySeqL = speD['Query_spe']
        #queryampliconL = []
        #for eachQueryMatch in querySeqL:
        #    id = eachQueryMatch[0]
        #    start,end = eachQueryMatch[2]
        #    queryampliconL.append(query_fastaD[id][start:end])
        #queryAmplicon = 'Query_spe'+':'+','.join(queryampliconL)
        #-----------------------------------------------------------
        #ampliconL = [queryAmplicon]
        #for otherspe in matchSpeL:
        #    targetSeqL = speD[otherspe]
        #    tmpL = []
        #    for eachQueryMatch in targetSeqL:
        #        id = eachQueryMatch[0]
        #        start,end = eachQueryMatch[2]
        #        tmpL.append(target_fastaD[id][start:end])
        #    targetAmplicon = otherspe+':'+','.join(tmpL)
        #    ampliconL.append(targetAmplicon)
            
        #----------------------------------------------------------------------

        outputL = [primer_id, '; '.join(matchSpeL), str(clone_spe)+'/'+str(targetSpeNum), str(differentAmpliconCnt), 
                ampliconsize, forward_primer, reverse_primer, ampliconFile]
        print >>finalOutput_fh, '\t'.join(outputL)
    
    #--------------------------------------------    
    finalOutput_fh.close()
#--------------------------------------------------------------------------

def readFasta(file, matchD={},gettargetSpeNum=0):
    seqD = {}
    speS = set()
    for line in open(file):
        if line[0] == '>':
            save = 0
            key = line.strip()[1:]
            if gettargetSpeNum:
                speS.add(key.split('___')[0])
            if matchD.get(key, ''):
                save = 1
                seqD[key] = []
        elif save:
            seqD[key].append(line.strip())
    #----------------------------------------
    for key,seqL in seqD.items():
        seqD[key] = ''.join(seqL)
    if gettargetSpeNum:
        return seqD, list(speS)
    return seqD
#------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    database = options.database
    ssr  = options.ssr
    flankMin, flankMax = [int(i) for i in options.flank.split(',')]
    am_size_opt, am_size_min, am_size_max = [int(i) for i in options.amplicon_size.split(',')]
    optsize, minsize, maxsize = [i for i in options.primer_size.split(',')]
    opttm, mintm, maxtm = [i for i in options.primer_tm.split(',')]
    output_prefix = options.op
    if not output_prefix:
        output_prefix = file
    log = output_prefix + '.log'
    log_fh = open(log, 'w')
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    misa_script = output_prefix + '.misa.pl'
    misa_ini    = 'misa.ini'
    if not ssr:
        outputMISA(misa_script)
        output_misa_ini(misa_ini)
        cmd = ' '.join(["perl", misa_script, file, misa_ini])
        if os.system(cmd):
            print >>sys.stderr, cmd + 'wrong'
            sys.exit(1)
        #--------------------------------
        ssr = file + '.misa'
    #------------------------------------------------------
    ssrDict = readmisaSSR(ssr)
    if ssrDict:
        print >>log_fh, "Get SSR"
    else:
        print >>log_fh, "No SSR found in given sequences. Please check your input and make sure the right sequences are given. If there is nothing wrong, mail to <y_yuan0732@163.com> for help."
        sys.exit(0)
    output_seq = output_prefix + '.ssr.flank.fa'
    extractFlankSeq(file, output_seq, ssrDict, flankMin, flankMax, am_size_min, am_size_max)
    
    #-----------------------------
    if (not os.path.exists(output_seq)) or os.stat(output_seq).st_size == 0:
        print >>log_fh, "No suitable SSR found due to too short flanking sequences. Please check your input and make sure the right sequences are given. If there is nothing wrong, mail to <y_yuan0732@163.com> for help."
        sys.exit(0)
    else:
        print >>log_fh, "Get flanking seq"
    #--------------------------------------
    all_primer_file = output_prefix + '.allPrimer'
    all_primer_file_fh = open(all_primer_file, 'w')
    deisgnPrimerForEachSeq(output_seq, all_primer_file_fh, optsize, minsize, maxsize, opttm, mintm, maxtm, 
        am_size_opt, am_size_min, am_size_max)
    all_primer_file_fh.close()
    
    #-----------------------------
    if (not os.path.exists(all_primer_file)) or os.stat(all_primer_file).st_size == 0:
        print >>log_fh, "No suitable primers found for amplicating SSR and their flanking regions. Please check your input and make sure the right sequences are given. If there is nothing wrong, mail to <y_yuan0732@163.com> for help."
        sys.exit(0)
    else:
        print >>log_fh, "Get flanking seq"
    #--------------------------------------

    primer_search_file     = ["primersearch -seqall", file, "-infile", all_primer_file,
            "-mismatchpercent 5", "-outfile", output_prefix+'.target.primerSearch']
    primer_search_file     = ' '.join(primer_search_file)
    if os.system(primer_search_file):
        print >>sys.stderr, primer_search_file + '**wrong'

    #---------------------------------------
    if (not os.path.exists(output_prefix+'.target.primerSearch')) or os.stat(output_prefix+'.target.primerSearch').st_size == 0:
        print >>log_fh, "Terrible wrong, please mail to <y_yuan0732@163.com> for help."
        sys.exit(0)
    else:
        print >>log_fh, "Get amplicon for search seq"
    #--------------------------------------

    primer_search_database = ["primersearch -seqall", database, "-infile", all_primer_file,
            "-mismatchpercent 5", "-outfile", output_prefix+'.database.primerSearch']
    primer_search_database     = ' '.join(primer_search_database)
    if os.system(primer_search_database):
        print >>sys.stderr, primer_search_database + '**wrong'

    #---------------------------------------
    if (not os.path.exists(output_prefix+'.database.primerSearch')) or os.stat(output_prefix+'.database.primerSearch').st_size == 0:
        print >>log_fh, "No suitable amplicons found in selected database. Double-check the right specie has been selected. If there is nothing wrong, mail to <y_yuan0732@163.com> for help."
        sys.exit(0)
    else:
        print >>log_fh, "Get amplicon for search seq"
    #--------------------------------------
    

    primerDict = {}
    querySeqD = readInPrimerSearch(output_prefix+'.target.primerSearch', primerDict, species="Query_spe")
    targetSeqD = readInPrimerSearch(output_prefix+'.database.primerSearch', primerDict)

    primerSeqD = dict([line.strip().split('\t',1) for line in open(all_primer_file)])
    query_fastaD = readFasta(file, targetSeqD)
    target_fastaD, targetSpeL = readFasta(database, targetSeqD, gettargetSpeNum=1)

    finalOutput = output_prefix+'.primer_evaluation.tsv'
    tabulizePrimerSearch(primerDict, primerSeqD, query_fastaD, target_fastaD, targetSpeL, finalOutput)
    

    #---------------------------------------
    if (not os.path.exists(finalOutput)) or os.stat(finalOutput).st_size == 0:
        print >>log_fh, "No primers have been found for discriminating input sequences among databases. If there is nothing wrong, mail to <y_yuan0732@163.com> for help."
        sys.exit(0)
    else:
        print >>log_fh, "Done"
    #--------------------------------------

    log_fh.close()
    """
    ssrDict = {'scaffold': 
                    [
                        {
                            ssr_type: "p1",
                            ssr: "(A)15",
                            size: 15,
                            start: 495,
                            end: 509
                        },
                        {
                            ssr_type: "p1",
                            ssr: "(A)15",
                            size: 15,
                            start: 495,
                            end: 509
                        },
                    ]
                }
    """




    #-------------------------------------------------------------------
#------------END main-----------------------------------

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " %         (' '.join(sys.argv), startTime, endTime)
    fh.close()
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------
#

