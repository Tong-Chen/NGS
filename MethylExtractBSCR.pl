#!/usr/bin/perl -w

###########################################################################
###########################################################################

#              *****************************************
#              *****  MethylExtract (version 1.8.1)  *****
#              ***************************************** 

#   		   Computational Genomics and Bioinformatics
#   		   Dept. of Genetics & Inst. of Biotechnology
#			  		   University of Granada

#              		Web: http://bioinfo2.ugr.es/

## This program is Copyright (C) 2012-13: Guillermo Barturen (bartg01@gmail.com), Antonio Rueda (aruemar@gmail.com), José L. Oliver (oliver@ugr.es), Michael Hackenberg (hackenberg@ugr.es)

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

# For questions, feedback, etc please contact to: Guillermo Barturen (bartg01@gmail.com), José L. Oliver (oliver@ugr.es, Michael Hackenberg (mlhack@gmail.com)

# To see the options, please launch MethylExtract without any command line arguments
############################################################################
############################################################################

use strict;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);

#################
#Default options#
#################
my %optionsDefault = (
	seqFile => "NA", inFile => "NA",
	FirstIgnor => 1000,
	LastIgnor => 1000,
	qscore => "phred33-quals", minQ => 20,
	flagW => "NA", flagC => "NA",
	tagW => "NA", tagC => "NA",
);
#################
#################

my %subQscore=("phred33-quals"=>\&sanger,"solexa-quals"=>\&solexa,"phred64-quals"=>\&illumina);

#Getting options
my ($inDir,$seqFile,$tagW,$tagC,$FirstIgnor,$qscore,$Q,$LastIgnor);
($inDir,$seqFile,$tagW,$tagC,$FirstIgnor,$qscore,$Q,$LastIgnor)=&GetOptions($tagW,$tagC,$FirstIgnor,$qscore,$Q,$LastIgnor);

print "\n############### Running MethylExtractBSCR v1.8.1 ###############\n\n";

&Main();

#############SUBPROCESS###################################################

sub GetOptions{
my %opts;
#Checking General Arguments & Help
if (@ARGV) {
	foreach (@ARGV) {
		my @GetOpt=split(/=/,$_);
		if ($GetOpt[0] eq "flagW" or $GetOpt[0] eq "flagC") {
			if ($GetOpt[1] eq "0") {$opts{$GetOpt[0]}="X"}
			else {$opts{$GetOpt[0]}=$GetOpt[1]}
		}
		elsif ($GetOpt[0] eq "tagW" or $GetOpt[0] eq "tagC") { 	# Deprecated
			if ($GetOpt[1] eq "0") {$opts{$GetOpt[0]}="X"}		# Deprecated
			else {$opts{$GetOpt[0]}=$GetOpt[1]}					# Deprecated
		}														# Deprecated
		else {$opts{$GetOpt[0]}=$GetOpt[1]}
	}
	#Checking input options
	foreach my $keyOpts (keys %opts){
		if (!$optionsDefault{$keyOpts}) {die "$keyOpts is not an accepted parameter, please check spelling and case sensitive\n";}
		else {}
	}
	#Checking sequences and indexed files paths
	if($opts{inFile}){
		$inDir=$opts{inFile};
		if (-e $inDir) {
			#Checking sam files
			my $checkAlign="N";
			if (-e $inDir) {
				$inDir =~ m/(\w+)$/;
				if ($1 eq "sam" or $1 eq "gz") {$checkAlign="Y";}
				else {}
			}
			else {}
			if ($checkAlign eq "N") {die "The input file doesn't seem to be *.sam or *.gz\n";}
			else {}
		}
		else {die "The files $inDir doesn't exist\n";}
	}
	else {die "Use inFile=[input file] to specify alignment file\n";}
	#Checking sequence path
	if($opts{seqFile}){
		$seqFile=$opts{seqFile};
		if (-e $seqFile) {}
		else {die "Cannot find Sequence File: $seqFile\n";}
	}
	else {die "Use seqFile=[sequence file] to set the fasta file to be used\n";}
	#Checking tags
	my (%tagsW,%tagsC);
	if ($opts{flagW}) {
		if ($opts{flagW} eq "X") {$tagsW{0}=1}
		else {%tagsW  = map { $_ => 1 } split(/,/, $opts{flagW})}
	}
	elsif ($opts{tagW}) {																	#	Deprecated
		if ($opts{tagW} eq "X") {$tagsW{0}=1}												#	Deprecated
		else {%tagsW  = map { $_ => 1 } split(/,/, $opts{tagW})}							#	Deprecated
		print "\nWarning: tagW is a deprecated parameter, it's been replaced by flagW\n";	#	Deprecated
	}																						#	Deprecated
	else {die "Watson FLAGs are required:\n Common single-end FLAG -> flagW=0\n Common pair-end FLAGs -> flagW=99,147\n"}
	if ($opts{flagC})  {
		if ($opts{flagC} eq "X") {$tagsC{0}=1}
		else {%tagsC = map { $_ => 1 } split(/,/, $opts{flagC})}
	}
	elsif ($opts{tagC}) {																	#	Deprecated
		if ($opts{tagC} eq "X") {$tagsW{0}=1}												#	Deprecated
		else {%tagsC  = map { $_ => 1 } split(/,/, $opts{tagC})}							#	Deprecated
		print "\nWarning: tagC is a deprecated parameter, it's been replaced by flagC\n";	#	Deprecated
	}																						#	Deprecated
	else {die "Crick FLAGs are required:\n Common single-end FLAG -> flagC=16\n Common pair-end FLAGs -> flagC=83,163\n"}
	$_[0]=\%tagsW;
	$_[1]=\%tagsC;
	my ($tagWout,$tagCout);
	foreach my $keys (keys %tagsW) {$tagWout.="$keys".","}
	chop($tagWout);
	foreach my $keys (keys %tagsC) {$tagCout.="$keys".","}
	chop($tagCout);
	print "Watson strand FLAG: $tagWout & Crick strand FLAG: $tagCout\n";
	#First number of positions ignored
	if (defined($opts{FirstIgnor})) {
		if ($opts{FirstIgnor}>=1000) {$_[2]=0}
		else {$_[2]=$opts{FirstIgnor}}
		if ($_[2]=~m/[0-9]+/) {
			$_[2]=~s/,/\./;
			if ($_[2]-int($_[2])>0) {die "Number of first bases ignored must be an integer\n";}
			else {}
		}
		else {die "Number of first bases ignored must be a real integer\n";}
		if ($_[2]==0) {print "First bases on the reads will not be ignored by default\n"}
		elsif ($_[2]>=10) {print "$_[2] first bases on reads to be ignored seems to be so high\n"}
		else {print "$_[2] first bases on reads will be ignored\n"}
	}
	else {
		print "First bases on the reads will not be ignored by default\n";
		$_[2]=0;
	}
	#Last number of positions ignored
	if (defined($opts{FirstIgnor})) {
		if ($opts{FirstIgnor}>=1000) {$_[5]=0}
		else {$_[5]=$opts{FirstIgnor}}
		if ($_[5]=~m/[0-9]+/) {
			$_[5]=~s/,/\./;
			if ($_[5]-int($_[5])>0) {die "Number of last bases ignored must be an integer\n";}
			else {}
		}
		else {die "Number of last bases ignored must be a real integer\n";}
		if ($_[5]==0) {print "Last bases on the reads will not be ignored by default\n"}
		elsif ($_[5]>=10) {print "$_[2] last bases on reads to be ignored seems to be so high\n"}
		else {print "$_[5] last bases on reads will be ignored\n"}
	}
	else {
		print "Last bases on the reads will not be ignored by default\n";
		$_[5]=0;
	}
	#Checking Qscore
	if ($opts{qscore}) {
		$qscore=$opts{qscore};
		if ($qscore eq "phred33-quals" or $qscore eq "phred64-quals" or $qscore eq "solexa-quals" or $qscore eq "NA"){}
		else {
			print "Quality score isn't an accepted format, launching by default phred33-quals\n";
			$_[3]="phred33-quals";
		}
	}
	else {$_[3]="phred33-quals";}
	#Minimun base quality
	if (defined($opts{minQ})) {
		$_[4]=$opts{minQ};
		if ($_[4]=~m/[0-9]+/) {
			$_[4]=~s/,/\./;
			if ($_[4]-int($_[4])>0) {die "Minimun base quality must be an integer\n";}
			elsif ($_[4]==0) {print "The base quality checker has been deactivated\n";}
			else {}
		}
		else {die "Minimun base quality must be an integer\n";}
	}
	else {$_[4]=20;}
}
else {
	print "\n################   MethylExtractBSCR   ###############\n";
	print "###############   Command-line help   ###############\n\n";
	print "Launch as:\n  perl MethExtractBSCR.pl seqFile=<sequence file> inFile=<alignments input file> flagW=<Watson FLAGs (multiple FLAGs comma separated)> flagC=<Crick FLAGs (multiple FLAGs comma separated)> [OPTIONS]\n";
	print "Optional Quality parameters:\n";
	print "  qscore=<fastq quality score: phred33-quals, phred64-quals,solexa-quals, solexa1.3-quals or NA> [default: phred33-quals]\n";
	print "  minQ=<minimun PHRED quality per sequenced nucleotide> [default: 20]\n  FirstIgnor=<number of first bases ignored> [default: 0]\n  LastIgnor=<number of last bases ignored> [default: 0]\n";
	die "\n";
}
return($inDir,$seqFile,$_[0],$_[1],$_[2],$_[3],$_[4],$_[5]);
}

sub Main {
	print "Reading $seqFile\n";
	my $seqst=&GetFasta("$seqFile");
	#extracting alignments
	my $countMeth=0;
	my $countAll=0;
	print "Checking bisulfite conversion rate\n";
	$inDir =~ m/(\w+)$/;
	if ($1 eq "sam") {
		open CHROM, "$inDir";
		while (my $line=<CHROM>) {
			if ($line!~/^@/i) {($countMeth,$countAll)=&conversionRate($line,$countMeth,$countAll,$seqst);}
			else {}
		}
		close CHROM;
	}
	elsif ($1 eq "gz") {
		my $z = new IO::Uncompress::AnyUncompress "$inDir";
		until (eof($z)) {
			my $line = <$z>;
			if ($line!~/^@/i) {($countMeth,$countAll)=&conversionRate($line,$countMeth,$countAll,$seqst);}
			else {}
		}
		close CHROM;
	}
	else {}
	my $bcr=1-($countMeth/$countAll);
	print "\n## MethylExtract_bisulfiteConversionRate Results ##\n";
	print "Number of cytosines: $countAll\n";
	print "Bisulfite conversion rate: $bcr\n"
}

#reading sequences
sub GetFasta {
	open (SEC, "$_[0]") or die "Fasta files do not seem to have chromosome alignments tags name\n";
  	my $seqst_temp = "";
    my $z = <SEC>;
    my $tes = substr($z,0,1);
    if($tes ne ">"){
     	die "Sequence seems not to be in fasta format !!!!";
    }
    $z=~ s/>//;
    $z =~ s/[\n\t\f\r\s]//g;
    while($z = <SEC>){
      	$z =~ s/[\n\t\f\r_0-9\s]//g;
      	$seqst_temp .= $z;
    }
    return($seqst_temp);
}

sub conversionRate {
	my @splitLine=split(/\t/,$_[0]);
	if ($tagW->{$splitLine[1]}) {
		#SingleEnd
		if ($splitLine[7]==0 and $splitLine[8]==0) {
			my $ini=$FirstIgnor;
			my $end=length($splitLine[9])-1-$LastIgnor;
			($_[1],$_[2])=&OrientationConversionRate(\@splitLine,$_[1],$_[2],$_[3],$ini,$end,"C","T");
		}
		#FirstPair
		elsif ($splitLine[7]>0 and $splitLine[8]>0) {
			my $ini=$FirstIgnor;
			my $end=length($splitLine[9])-1-$LastIgnor;
			($_[1],$_[2])=&OrientationConversionRate(\@splitLine,$_[1],$_[2],$_[3],$ini,$end,"C","T");
		}
		#SecondPair
		elsif ($splitLine[7]>0 and $splitLine[8]<0) {
			my $ini=$LastIgnor;
			my $end=length($splitLine[9])-1-$FirstIgnor;
			($_[1],$_[2])=&OrientationConversionRate(\@splitLine,$_[1],$_[2],$_[3],$ini,$end,"C","T");
		}
		else {}
	}
	elsif ($tagC->{$splitLine[1]}) {
		#SingleEnd
		if ($splitLine[7]==0 and $splitLine[8]==0) {
			my $ini=$LastIgnor;
			my $end=length($splitLine[9])-1-$FirstIgnor;
			($_[1],$_[2])=&OrientationConversionRate(\@splitLine,$_[1],$_[2],$_[3],$ini,$end,"G","A");
		}
		#FirstPair
		elsif ($splitLine[7]>0 and $splitLine[8]>0) {
			my $ini=$FirstIgnor;
			my $end=length($splitLine[9])-1-$LastIgnor;
			($_[1],$_[2])=&OrientationConversionRate(\@splitLine,$_[1],$_[2],$_[3],$ini,$end,"G","A");
		}
		#SecondPair
		elsif ($splitLine[7]>0 and $splitLine[8]<0) {
			my $ini=$LastIgnor;
			my $end=length($splitLine[9])-1-$FirstIgnor;
			($_[1],$_[2])=&OrientationConversionRate(\@splitLine,$_[1],$_[2],$_[3],$ini,$end,"G","A");
		}
		else {}
	}
	else {}
	return($_[1],$_[2]);
}

sub OrientationConversionRate {
	for (my $i=$_[4];$i<=$_[5];$i++) {
			if (uc(substr($_[3],${$_[0]}[3]+$i-1,1)) eq "$_[6]") {
				my $ascii=substr(${$_[0]}[10],$i,1);
				my $PHREDval=&{$subQscore{$qscore}}($ascii);
				if ($PHREDval>=$Q) {
					if (uc(substr(${$_[0]}[9],$i,1)) eq "$_[6]") {
						$_[1]++;
						$_[2]++;
					}
					elsif (uc(substr(${$_[0]}[9],$i,1)) eq "$_[7]") {$_[2]++}
					else {}
				}
				else {}				
			}
			else {}
	}
	return($_[1],$_[2]);
}

sub QScore {
	my $ReturnData;
	if ($_[0]>=33 and $_[0]<59) {
		#$ReturnData="fastq-sanger";
		$ReturnData="phred33-quals";
	}
	elsif ($_[0]>=59 and $_[0]<64) {
		#$ReturnData="fastq-solexa";
		$ReturnData="solexa-quals";
	}
	elsif ($_[0]>=64 and $_[0]<=126) {
		#$ReturnData="fastq-illumina";
		$ReturnData="phred64-quals";
	}
	else {die "Cannot stablish Quality Score format\n"}
	return($ReturnData);
}

#Qscore calculate
sub sanger {
	return(ord($_[0])-33);
}
sub solexa {
	return((10 * log(1 + 10 ** (ord($_[0]) - 64) / 10)) / log(10));
}
sub illumina {
	return(ord($_[0])-64);
}
