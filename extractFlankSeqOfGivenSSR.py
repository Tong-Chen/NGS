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
    Extract flanking sequences of given regions.

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
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

def outputMISA(misa="misa.pl"):
    if os.path.exists(misa):
        return

    misa_fh = open(misa,w)
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


#§§§§§ DECLARATION §§§§§#

# Check for arguments. If none display syntax #

if (@ARGV == 0)
  {
  open (IN,"<$0");
  while (<IN>) {if (/^\#\# (.*)/) {$message .= "$1\n"}};
  close (IN);
  die $message;
  };

# Check if help is required #

if ($ARGV[0] =~ /-help/i)
  {
  open (IN,"<$0");
  while (<IN>) {if (/^\#\#\#(.*)/) {$message .= "$1\n"}};
  close (IN);
  die $message;
  };

# Open FASTA file #

open (IN,"<$ARGV[0]") || die ("\nError: FASTA file doesn't exist !\n\n");
open (OUT,">$ARGV[0].misa");
print OUT "ID\tSSR nr.\tSSR type\tSSR\tsize\tstart\tend\n";

# Reading arguments #

open (SPECS,"misa.ini") || die ("\nError: Specifications file doesn't exist !\n\n");
my %typrep;
my $amb = 0;
while (<SPECS>)
   {
   %typrep = $1 =~ /(\d+)/gi if (/^def\S*\s+(.*)/i);
   if (/^int\S*\s+(\d+)/i) {$amb = $1}
   };
my @typ = sort { $a <=> $b } keys %typrep;


#§§§§§ CORE §§§§§#

$/ = ">";
my $max_repeats = 1; #count repeats
my $min_repeats = 1000; #count repeats
my (%count_motif,%count_class); #count
my ($number_sequences,$size_sequences,%ssr_containing_seqs); #stores number and size of all sequences examined
my $ssr_in_compound = 0;
my ($id,$seq);
while (<IN>)
  {
  next unless (($id,$seq) = /(.*?)\n(.*)/s);
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
    my $search = "(([acgt]{$motiflen})\\2{$minreps,})";
    while ( $seq =~ /$search/ig ) #scan whole sequence for that class
      {
      my $motif = uc $2;
      my $redundant; #reject false type motifs [e.g. (TT)6 or (ACAC)5]
      for ($j = $motiflen - 1; $j > 0; $j--)
        {
        my $redmotif = "([ACGT]{$j})\\1{".($motiflen/$j-1)."}";
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
    print OUT "$id\t$count_seq\t$ssrtype\t$ssrseq\t",($end - $start + 1),"\t$start\t$end\n"
    };
  };

close (OUT);
open (OUT,">$ARGV[0].statistics");

#§§§§§ INFO §§§§§#

#§§§ Specifications §§§#
print OUT "Specifications\n==============\n\nSequence source file: \"$ARGV[0]\"\n\nDefinement of microsatellites (unit size / minimum number of repeats):\n";
for ($i = 0; $i < scalar (@typ); $i++) {print OUT "($typ[$i]/$typrep{$typ[$i]}) "};print OUT "\n";
if ($amb > 0) {print OUT "\nMaximal number of bases interrupting 2 SSRs in a compound microsatellite:  $amb\n"};
print OUT "\n\n\n";

#§§§ OCCURRENCE OF SSRs §§§#

#small calculations
my @ssr_containing_seqs = values %ssr_containing_seqs;
my $ssr_containing_seqs = 0;
for ($i = 0; $i < scalar (@ssr_containing_seqs); $i++) {$ssr_containing_seqs += $ssr_containing_seqs[$i]};
my @count_motifs = sort {length ($a) <=> length ($b) || $a cmp $b } keys %count_motifs;
my @count_class = sort { $a <=> $b } keys %count_class;
for ($i = 0; $i < scalar (@count_class); $i++) {$total += $count_class{$count_class[$i]}};

#§§§ Overview §§§#
print OUT "RESULTS OF MICROSATELLITE SEARCH\n================================\n\n";
print OUT "Total number of sequences examined:              $number_sequences\n";
print OUT "Total size of examined sequences (bp):           $size_sequences\n";
print OUT "Total number of identified SSRs:                 $total\n";
print OUT "Number of SSR containing sequences:              $ssr_containing_seqs\n";
print OUT "Number of sequences containing more than 1 SSR:  ",$ssr_containing_seqs - ($ssr_containing_seqs{1} || 0),"\n";
print OUT "Number of SSRs present in compound formation:    $ssr_in_compound\n\n\n";

#§§§ Frequency of SSR classes §§§#
print OUT "Distribution to different repeat type classes\n---------------------------------------------\n\n";
print OUT "Unit size\tNumber of SSRs\n";
my $total = undef;
for ($i = 0; $i < scalar (@count_class); $i++) {print OUT "$count_class[$i]\t$count_class{$count_class[$i]}\n"};
print OUT "\n";

#§§§ Frequency of SSRs: per motif and number of repeats §§§#
print OUT "Frequency of identified SSR motifs\n----------------------------------\n\nRepeats";
for ($i = $min_repeats;$i <= $max_repeats; $i++) {print OUT "\t$i"};
print OUT "\ttotal\n";
for ($i = 0; $i < scalar (@count_motifs); $i++)
  {
  my $typ = length ($count_motifs[$i]);
  print OUT $count_motifs[$i];
  for ($j = $min_repeats; $j <= $max_repeats; $j++)
    {
    if ($j < $typrep{$typ}) {print OUT "\t-";next};
    if ($count_motifs[$i]->{$j}) {print OUT "\t$count_motifs[$i]->{$j}"} else {print OUT "\t"};
    };
  print OUT "\t$count_motifs{$count_motifs[$i]}\n";
  };
print OUT "\n";

#§§§ Frequency of SSRs: summarizing redundant and reverse motifs §§§#
# Eliminates %count_motifs !
print OUT "Frequency of classified repeat types (considering sequence complementary)\n-------------------------------------------------------------------------\n\nRepeats";
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
for ($i = $min_repeats; $i <= $max_repeats; $i++) {print OUT "\t$i"};
print OUT "\ttotal\n";
@red_rev = sort {length ($a) <=> length ($b) || $a cmp $b } keys %red_rev;
for ($i = 0; $i < scalar (@red_rev); $i++)
  {
  my $typ = (length ($red_rev[$i])-1)/2;
  print OUT $red_rev[$i];
  for ($j = $min_repeats; $j <= $max_repeats; $j++)
    {
    if ($j < $typrep{$typ}) {print OUT "\t-";next};
    if ($red_rev[$i]->{$j}) {print OUT "\t",$red_rev[$i]->{$j}}
    else {print OUT "\t"}
    };
  print OUT "\t",$red_rev[$i]->{"total"},"\n";
  };

"""
    misa_fh.close()
#----outputMISA---------------------------------------------------------------------

def output_misa_ini(misa_ini="misa.ini"):
    if os.path.exists(misa_ini):
        return

    misa_ini_fh = open(misa_ini,w)
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
    parser.add_option("-s", "--ssr-file", dest="ssr",
        help="SSR file generated by misa.pl. If not given, SSR will be first identified.")
    parser.add_option("-S", "--amplicon-size", dest="amplicon_size",
        default = "100,280", 
        help="Amplicon size. Default <100,280> meaning at least 100 nt and at most 280 nt will be extracted.")
    parser.add_option("-f", "--flank", dest="flank",
        default = "13,250", 
        help="Length of flank regions to be extracted along each SSR. Default <13,250> meaning at least 13 nt and at most 250 nt will be extracted.")
    parser.add_option("-p", "--output-prefix", dest="op",
        help="Output prefix.")
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
        print >>output_seqfh, ">{}\t{}\t{}\t{}".format(key, ssr_type, ssr, size)
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
        if line[0]== '>':
            if seqL:
                seq = ''.join(seqL)
                ssrL = ssrDict.get(key,[])
                if not ssrL:
                    key2 = key.replace(' ','_')
                    ssrL = ssrDict.get(key2,[])
                if ssrL:
                    extractSeq(output_seqfh,key,seq,ssrL,flankMin,flankMax, am_size_min, am_size_max)
            key = line[1:-1]
            seqL = []
        else:
            seqL.append(line.strip())
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

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    ssr  = options.ssr
    flankMin, flankMax = [int(i) for i in options.flank.split(',')]
    am_size_min, am_size_max = [int(i) for i in options.amplicon_size.split(',')]
    output_prefix = options.op
    if not output_prefix:
        output_prefix = file
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    if not ssr:
        outputMISA()
        output_misa_ini()
        cmd = ' '.join(["perl misa.pl ", file, "misa.ini"])
        if os.system(cmd):
            print >>sys.stderr, cmd + 'wrong'
            sys.exit(1)
        #--------------------------------
        ssr = file + '.misa'
    #------------------------------------------------------
    ssrDict = readmisaSSR(ssr)

    output_seq = output_prefix + '.ssr.flank.fa'
    extractFlankSeq(file, output_seq, ssrDict, flankMin, flankMax, am_size_min, am_size_max)

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

