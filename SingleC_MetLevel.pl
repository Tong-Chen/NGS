#!/usr/bin/perl -w
use strict;

if(@ARGV != 2){
	print "Usage: perl $0 <reference> <pileup or - represents stdin> > <out>\n";
	exit;
}

my $ref = $ARGV[0];
my %ref;
open IN,"$ref" or die $!;#genome.fa
my ($Chr,$seq) = ('','');
while(<IN>)
{
        chomp;
        if (/>/){
                $_ =~ />(.*)/;
                my $tmp_chr = $1;
                if ($tmp_chr ne $Chr){
                        if ($seq ne ''){
                                $ref{$Chr} = $seq;
                                $seq = '';
                        }
                        $Chr = $tmp_chr;
                }
                next;
        }
        $seq .= $_;
}
$ref{$Chr} = $seq;
close IN;

my $in;
print "#Chr\tPos\tRef\tChain\tTotal\tMet\tUnMet\tMetRate\tRef_context\tType\n";
my %hash;
#if ("$ARGV[1]" == "-") {
#	$in = *STDIN;
#} else {
open IN,"$ARGV[1]" or die $!;
#}
while(<IN>){
	chomp;
	my @a = split /\s+/,$_;
	my $chr = $a[0];
	my $pos = $a[1];
	my $ref = $a[2];
	my $depth = $a[3];
	next if ($ref !~ /[CG]/i or $chr =~ /random/ or $chr eq "chrM" or $depth == 0);
	my $tmp_pos = $pos - 1;
	my ($total,$rate,$chain,$dot,$base,$str,$Ctype);
	if ($ref =~ /C/i){
		$dot = $a[4] =~ tr/\./\./;
		$base = $a[4] =~ tr/T/T/;
		$total = $dot + $base;
		next if ($total == 0);
		$chain = "+";
		$str = substr($ref{$chr},$tmp_pos,3);
		$str = uc $str;
		if (length($str) == 2){
			my $base2 = substr($str,1,1);
			if ($base2 =~ /G/i){
                                $Ctype = "CpG";
                        }
		}
		else{
                        if($str=~/C[ATC]G/i){
				$Ctype = "CHG";
                        }
                        elsif ($str =~ /C[ATC][ATC]/){
                                $Ctype = "CHH";
                        }
                        elsif ($str =~ /^CG/i){
                                $Ctype = "CpG";
                        }
		}
	}
	elsif ($ref =~ /G/i){
		$dot = $a[4] =~ tr/,/,/;
		$base = $a[4] =~ tr/a/a/;
		$total = $dot + $base;
		next if ($total == 0);
		$chain = "-";
                next if ($tmp_pos == 0);
		if ($tmp_pos == 1){
                        $str = substr($ref{$chr},0,2);
                }
                elsif ($tmp_pos > 1){
                        $tmp_pos = $tmp_pos - 2;
                        $str = substr($ref{$chr},$tmp_pos,3);
                }
                $str = reverse $str;
                $str = uc $str;
                $str =~ tr/ATGC/TACG/;
                if (length($str) == 2){
                        my $base2 = substr($str,1,1);
                        if ($base2 =~ /G/i){
                                $Ctype = "CpG";
                        }
                }
                else{
                        if($str=~/C[ATC]G/i){
                                $Ctype = "CHG";
                        }
                        elsif ($str =~ /C[ATC][ATC]/){
                                $Ctype = "CHH";
                        }
                        elsif ($str =~ /^CG/i){
                                $Ctype = "CpG";
                        }
                }
	}
	$rate = $dot/$total;
	if (defined $Ctype){
		print "$chr\t$pos\t$ref\t$chain\t$total\t$dot\t$base\t$rate\t$str\t$Ctype\n";
	}
}

#if ("$ARGV[1]" != "-") {
#close $in;
#}
close IN;
