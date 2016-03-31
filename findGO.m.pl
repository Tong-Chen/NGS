#!/usr/bin/env perl
use warnings;
use lib "/MPATHB/soft/homer/.//bin";
my $homeDir = "/MPATHB/soft/homer/./";
my $promoterSeqOffset = -2000;

my $goDir = $homeDir . "/data/GO/";
my $goFile = $homeDir . "/data/GO/GO.txt";
my $promoterDir = $homeDir . "/data/promoters/";
my $accDir = $homeDir . "/data/accession/";
my $pvalueThresh = 0.05;


sub printCMD {
	print STDERR "\n\tUsage: findGO.pl <target ids file> <organism> <output directory> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-bg <background ids file>\n";
	print STDERR "\t\t-cpu <#> (number of cpus to use)\n";
	print STDERR "\t\t-human (convert IDs and run as human [uses homologene])\n";
	print STDERR "\n";
	exit;
}
if (@ARGV < 3) {
	printCMD();
}
my $inputIDfile = $ARGV[0];
my $org = $ARGV[1];
my $outDir = $ARGV[2];
my $bgIDfile = '';
my $maxCPUs = 1;
my $humanFlag = 0;

for (my $i=3;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-bg') {
		$bgIDfile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-p' || $ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-human') {
		$humanFlag = 1;
	} else {
		printCMD();
	}
}

$tmpFileFg = rand() . ".fg.tmp";
$tmpFileBg = rand() . ".bg.tmp";

if ($humanFlag == 0) {
	`convertIDs.pl "$inputIDfile" $org gene > "$tmpFileFg"`;
	if ($bgIDfile ne '') {
		`convertIDs.pl "$bgIDfile" $org gene > "$tmpFileBg"`;
	} else {
		if (-f "$promoterDir/$org.base.gene") {
			`cp "$promoterDir/$org.base.gene" "$tmpFileBg"`;
		}
	}
} else {
	`convertOrganismID.pl "$inputIDfile" $org human gene > "$tmpFileFg"`;
	if ($bgIDfile ne '') {
		`convertOrganismID.pl "$bgIDfile" $org human gene > "$tmpFileBg"`;
	} else {
		if (-f "$promoterDir/human.base.gene") {
			`cp "$promoterDir/human.base.gene" "$tmpFileBg"`;
		}
	}
	$org = 'human';
}

%geneNames = ();
open IN, $accDir . $org . ".description";
while (<IN>) {
	chomp;
	my @line = split /\t/;
	my $name = $line[4];
	if (@line < 5 || $line[4] eq '-' || $line[4] eq '') {
		$name = $line[0];
	} 
	$geneNames{$line[0]} = $name;
}
close IN;

my @ontologies = ();
open IN, $goFile;
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^#/);
	my @line = split /\t/;
	my $go = {name=>$line[0],genes=>$line[1],output=>$line[2],desc=>$line[3],url=>$line[4],linkName=>$line[5]};
	push(@ontologies,$go);
	
}
close IN;


`mkdir -p "$outDir"`;

my $cpus=0;
foreach(@ontologies) {
	my $go = $_;
	my $pid = fork();
	$cpus++;
	if ($pid == 0) {
		#child process
		`findGOtxt.pl "$goDir/$org.$go->{'genes'}" "$tmpFileBg" "$tmpFileFg" > "$outDir/$go->{'output'}"`;
		exit(0);
	}
	if ($cpus >= $maxCPUs) {
		wait();
		$cpus--;
	}
}
my $id = 0;
while ($id >= 0) {
	$id = wait();
}
`/bin/rm -f "$tmpFileFg" "$tmpFileBg"`;	

my %results = ();
for (my $i=0;$i<@ontologies;$i++) {
	my $go = $ontologies[$i];
	parseResults("$outDir/$go->{'output'}", $go->{'name'}, $pvalueThresh);
}

#my @go = sort {$results{$a}->{'lp'} <=> $results{$b}->{'lp'}} keys %results;
#
#
#open MAIN, ">$outDir/geneOntology.html";
#print MAIN "<HTML><HEAD><TITLE>Gene Ontology Results</TITLE></HEAD><BODY>\n";
#print MAIN "<H1>Gene Ontology Enrichment Results</H1>\n";
#print MAIN "<H4>Text file version of complete results (i.e. open with Excel)\n";
#print MAIN "<UL>\n";
#
#foreach(@ontologies) {
#	my $go = $_;
#	
#	print MAIN "<LI><A HREF=\"$go->{'output'}\">$go->{'name'}</A>: $go->{'desc'} ";
#	print MAIN "(<A HREF=\"$go->{'url'}\">$go->{'linkName'}</A>)</LI>\n";
#
#}
#
#print MAIN "</UL>\n";
#print MAIN "<H2>Enriched Categories</H2>\n";
#print MAIN "<TABLE border=\"1\" cellpading=\"0\" cellspacing=\"0\">\n";
#print MAIN "<TR><TH>P-value</TH><TD>ln(P)</TD><TD>Term</TD><TD>GO Tree</TD><TD>GO ID</TD>";
#print MAIN "<TD># of Genes in Term</TD><TD># of Target Genes in Term</TD><TD># of Total Genes</TD>";
#print MAIN "<TD># of Target Genes</TD><TD>Common Genes</TD></TR>\n";
#
#foreach(@go) {
#	my $goid = $_;
#	my $pvalue = sprintf("%.3e", $results{$goid}->{'p'});
#	last if ($pvalue > $pvalueThresh);
#
#	my $lp = sprintf("%.2f", $results{$goid}->{'lp'});
#	my $term = $results{$goid}->{'term'};
#	my $tree = $results{$goid}->{'n'};
#	my $nterm = $results{$goid}->{'t'};
#	my $ncommon = $results{$goid}->{'c'};
#	my $N = $results{$goid}->{'N'};
#	my $T = $results{$goid}->{'P'};
#
#	print MAIN "<TR><TD>$pvalue</TD><TD>$lp</TD><TD>$term</TD><TD>$tree</TD><TD>$goid</TD>";
#	print MAIN "<TD>$nterm</TD><TD>$ncommon</TD><TD>$N</TD><TD>$T</TD>\n";
#	print MAIN "<TD>";
#	my $z = 0;
#	foreach(@{$results{$goid}->{'genes'}}) {
#		print MAIN "," if ($z > 0);
#		$z++;
#		my $name = $_;
#		#if (exists($geneNames{$_})) {
#		#	$name = $geneNames{$_};
#		#}
#		print MAIN "$name";
#	}
#	print MAIN "</TD></TR>\n";
#
#}
#print MAIN "</TABLE>\n";
#print MAIN "</BODY></HTML>\n";
#close MAIN;


sub parseResults {
	my ($file, $name, $pvalueThresh) = @_;
	my $c = 0;
	my $output = 0;
	open IN, $file;
	while (<IN>) {
		$c++;
		next if ($c < 2);
		chomp;
		my @line = split /\t/;


		my $term = $line[1];
		my $goid = $line[0];
		my $pvalue = $line[2];
		my $fdr = $pvalue;
		my $lp = $line[3];
		my $numTerm = $line[4];
		my $common = $line[5];
		my $total = $line[8];
		my $pos = $line[7];
		my @common = ();
		if (@line > 9) {
			my @gids = split /\,/, $line[9];
			foreach(@gids) {
				if (exists($geneNames{$_})) {
					push(@common, $geneNames{$_});
				} else {
					push(@common, $_);
				}
			}
		}
		#$results{$goid} = {term=>$term, lp=>$lp, p=>$pvalue, 
		#	fdr=>$fdr, t=>$numTerm,c=>$common,
		#	N=>$total,P=>$pos,n=>$name,genes=>\@common};

		if ($pvalue <= $pvalueThresh){
			if ($output == 0 ){
				open TMP, ">$tmpFileFg";
				print TMP "TermID\tTerm\tFDR\tneg_log10FDR\tGenes in Term\tTarget Genes in Term\tFraction of Targets in Term\tTotal Target Genes\tTotal Genes\tGene Symbols\tEntrez Gene IDs\n";
				$output = 1;
			}

			print TMP "$line[0]";
			#for (my $i=1;$i<3;$i++) {
				#print TMP "\t$line[$i]";
			#}
			#print TMP "\t$pvalue";
			for (my $i=1;$i<9;$i++) {
				print TMP "\t$line[$i]";
			}
			if (scalar(@common) > 0) {
				print TMP "\t$common[0]";
				for (my $i=1;$i<@common;$i++) {
					print TMP ",$common[$i]";
				}
			} else {
				print TMP "\tNULL"
			}
			for (my $i=9;$i<@line;$i++) {
				print TMP "\t$line[$i]";
			}
			print TMP "\n";
		}
	}
	close IN;
	if ($output == 1){
		close TMP;
		`multipleTestCorrection.sh -f $tmpFileFg -k FALSE -M FDR -L neg_log10FDR -o $file`
	} else {
		`/bin/rm -f "$file"`
	}
	
	#`mv "$tmpFileFg" "$file"`;
}


