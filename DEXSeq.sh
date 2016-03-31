#!/bin/bash

#set -x
set -e
set -u

usage()
{
cat <<EOF >&2
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to perform DE exon analysis using DEXSeq.

It requires at least two input files.

The count file (normally generated using ${txtred}HTseq-count${txtrst} and multiple
samples are pasted using ${txtred}pasteMultipleFilesSpecialCol.py${txtrst}):

#--------FILE content-----------------------------
ID	A_1	A_2	A_3	B_1	B_2	C_1	C_2	...
a	1	1	1	2	2	3	3	...
b	1	1	1	2	2	3	3	...
c	1	1	1	2	2	3	3	...
d	1	1	1	2	2	3	3	...
e	1	1	1	2	2	3	3	...
#--------FILE content-----------------------------

The sample file (with the first column as the colnames of countfile
and second column indicates the origin of each replicates.)

#--------FILE content------------------------------
#---the formula should be <conditions> which is the default
Sample	conditions
A_1	A
A_2	A
A_3	A
B_1	B
B_2	B
C_1	C
C_2	C
#--------FILE content-----------------------------
#--------FILE content--multiple columns are allowed-------------
#--------One <conditions> column is needed----------------------
#---the formula should be <cell+time+cell:time> 
#---the reduced formula (only applicable in <time-series>) 
#---should be  <cell+time> or <cell> or <time> for
#---specific purposes
#---First appeared sample will be treated as control.
Sample	conditions	cell	time
A_1_1	A_1	A	1	
A_1_2	A_1	A	1
A_1_3	A_1	A	1
A_2_1	A_2	A	2	
A_2_2	A_2	A	2
A_2_3	A_2	A	2
B_1_1	B_1	B	1	
B_1_2	B_1	B	1
B_1_3	B_1	B	1
B_2_1	B_2	B	2	
B_2_2	B_2	B	2
B_2_3	B_2	B	2
#--------FILE content-----------------------------

Optional file to specify the comparasion patterns you want to do


#--------FILE content- for pairwise compare-------
A	B
A	C
B	C
#--------FILE content-----------------------------

#--------FILE content- for one control-------
A	C
B	C
#--------FILE content-----------------------------

Entrez file
#--------FILE content-----------------------------
ID	EntrezGene
AT3G18710	821402
AT4G25880	828694
#--------FILE content-----------------------------


${txtbld}OPTIONS${txtrst}:
	-f	Directory containing output of <dexseq_count.py>. 
		The suffix for count files should be <.DEXSeq>.
		${bldred}[A directory containing count files ended with 
		<.DEXSeq>,  NECESSARY]
		${txtrst}
	-Q	The prefix of for output files.
		${bldred}NECESSARY${txtrst}
	-s	Sample file ${bldred}[A multiple columns file with header line, 
		For <timeseries>,  one <conditions> columns is needed.
		NECESSARY]
		CHECK ABOVE FOR DETAILS
		${txtrst}
	-g	The GTF file used for <dexseq_count.py>.
	-d	The design formula for DEXSeqDataSetFromMatrix.
		${bldred}[Default <conditions>, 
		accept <cell+time+cell:time> for example 2.]
		${txtrst}
	-D	The reduced design formula for DEXSeq.
		${bldred}[Only applicable to <timeseries> analysis, 
		accept <cell+time> or <time> or <cell> for example 2.]
		${txtrst}
	-a	All sample mix analysis.
		${bldred}[Default <TRUE>, accept <FALSE>.${txtrst}
	-m	Specify the comparasion mode.
		${bldred}[Default <pairwise>, accept <timeseries>,
		<pairwise> comparasion will still be done in <timeseries>
		mode.
		NECESSARY]
		${txtrst}
	-p	A file containing the pairs needed to do comparasion. 
		CHECK ABOVE FOR DETAILS
		All samples will be compared in <pairwise> mode if not specified here.
	
	-F	Log2 Fold change for screening DE genes.
		${bldred}Default 1${txtrst}

	-P	FDR for screening DE genes.
		${bldred}Default 0.01${txtrst}

	-q	FDR for screening time-series DE genes.
		${bldred}Default 0.1${txtrst}

	-e	Execute programs
		${bldred}Default TRUE${txtrst}

	-i	Install packages of not exist. 
		${bldred}Default FALSE${txtrst}

Eg.
	$0 -f "'dir/1.txt','dir/2.txt'" -Q dir/prefix.dir -s sampleFile -g gtf -d conditions
	$0 -f "'dir/1.txt','dir/2.txt'" -Q dir/prefix.dir -s sampleFile -g gtf -p compare_pair
	$0 -f "'dir/1.txt','dir/2.txt'" -Q dir/prefix.dir -s sampleFile -g gtf -m timeseries -d cell+time+cell:time -D cell+time -p compare_pair

EOF
}

dir=
prefix=
sample=
formula="conditions"
reducedFormula="nouse"
all='TRUE'
compare_mode='pairwise'
compare_pair='FALSE'
header='TRUE'
ist='FALSE'
execute='TRUE'
#outputdir='./'
fdr=0.01
log2fc=1
t_fdr=0.1
gtf=

while getopts "ha:d:D:e:f:F:g:i:m:o:p:P:q:Q:s:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		a)
			all=$OPTARG
			;;
		d)
			formula=$OPTARG
			;;
		D)
			reducedFormula=$OPTARG
			;;
		e)
			execute=$OPTARG
			;;
		f)
			dir=$OPTARG
			;;
		F)
			log2fc=$OPTARG
			;;
		g)
			gtf=$OPTARG
			;;
		i)
			ist=$OPTARG
			;;
		m)
			compare_mode=$OPTARG
			;;
		o)
			outputdir=$OPTARG
			;;
		p)
			compare_pair=$OPTARG
			;;
		P)
			fdr=$OPTARG
			;;
		q)
			t_fdr=$OPTARG
			;;
		Q)
			prefix=$OPTARG
			;;
		s)
			sample=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

if [ -z "$dir" ]; then
	usage
	exit 1
fi

mid=".DEXSeq"

formulaV=`echo ${formula} | sed 's/ //g' | awk '{a=split($1,b,"+"); c="\""b[1]"\""; for (i=2;i<=a;i++) if (b[i] !~/:/) c=c",""\""b[i]"\""; print c;}'`
first_v=`echo ${formula} | sed 's/ //g' | cut -d ':' -f 1 | cut -d '+' -f 1`
second_v=`echo ${formula} | sed 's/ //g' | cut -d ':' -f 1 | cut -d '+' -f 2`

#formulaS=`echo ${formula} | sed 's/ //g' | awk '{a=split($1,b,"+"); c=b[1]; for (i=2;i<=a;i++) c=c":"b[i]; print $1"+"c;}'`

all_de="${prefix}${mid}.all.DE"

cat <<END >${prefix}${mid}.r

if ($ist){
	source("https://bioconductor.org/biocLite.R")
	source(pipe(paste("wget -O -", URL)))
	biocLite("DEXSeq")
	biocLite("BiocParallel")
}

library(DEXSeq)
library("RColorBrewer")
library("gplots")
library("amap")
library("ggplot2")
library("BiocParallel")

#data <- read.table("${dir}/${prefix}", header=T, row.names=1, com='', quote='',
#	check.names=F, sep="\t")

#data <- data[rowSums(data)>2,]

#countFiles <- list.files("${dir}", pattern=".DEXSeq\$", full.names=T)
gffFile <- "${gtf}"

sample <- read.table("${sample}", header=T, row.names=1, com='',
	quote='', check.names=F, sep="\t", colClasses="factor")

sample_rowname <- rownames(sample)

sample <- data.frame(lapply(sample, function(x) factor(x, levels=unique(x))))
rownames(sample) <- sample_rowname


#paste0(sample \$cell, sample\$time)

#if ("${compare_mode}" == "pairwise") {
#	print("Perform pairwise comparasion using <design=~${formula}>")
#	ddsFullCountTable <- DEXSeqDataSetFromHTSeq(countFiles,
#		sampleData = sample,  
#		design= ~ sample + exon + ${formula}:exon,
#	   	flattenedfile=gffFile)
#} else if ("${compare_mode}" == "timeseries") {
#	print("Perform pairwise comparasion using <design=~conditions>")
#	ddsFullCountTable <- DEXSeqDataSetFromMatrix(countData = data,
#		colData = sample,  design= ~conditions)
#}

DEXSeq_ct <- function
(sample, gffFile, prefix, comp, fdr
){
	sample_rowname <- rownames(sample)
	countFiles <- file.path("${dir}", paste0(sample_rowname, ".DEXSeq"))
	ddsFullCountTable <- DEXSeqDataSetFromHTSeq(countFiles,
		sampleData = sample,  
		design= ~ sample + exon + ${formula}:exon,
		flattenedfile=gffFile)

	BPPARAM <- BiocParallel::MulticoreParam(workers=30)

	print("Normalization")
	dxd <- estimateSizeFactors(ddsFullCountTable)

	print("Dispersion estimation")
	dxd <- estimateDispersions(dxd, BPPARAM=BPPARAM)
	print("Testing for differential exon usage")
	dxd <- testForDEU(dxd, BPPARAM=BPPARAM)
	dxd <- estimateExonFoldChanges(dxd, fitExpToVar="${formula}",
		BPPARAM=BPPARAM)
	dds <- DEXSeqResults(dxd)

	# Get normalized counts
	print("Output normalized counts")
	normalized_counts <- counts(dds, normalized=TRUE)
	
	file = paste0(prefix, ".normalized.xls")
	write.table(normalized_counts, file=file,
		quote=F, sep="\t", row.names=T, col.names=T)
	system(paste("sed -i '1 s/^/ID\t/'", file))

	DEXSeqHTML(dds,path=prefix, FDR=fdr, fitExpToVar="${formula}",
	BPPARAM=BPPARAM, file=paste0(comp, ".DEU.html"))

	file = paste0(prefix, ".results.xls")
	write.table(dds, file=file, quote=F, sep="\t", row.names=T,
		col.names=T)	
	system(paste("sed -i '1 s/^/ID\t/'", file))
}

de_twosample <- function
(
sample,
gffFile,
fdr,
sampleV
){
	sampA <- as.vector(sampleV\$sampA)
	sampB <- as.vector(sampleV\$sampB)
	comp <- paste(sampA, "_vs_", sampB, sep=".")
	prefix_name <- paste("${prefix}${mid}",comp,sep=".") 
	sampleTmp<-sample[sample\$conditions==sampA|sample\$conditions==sampB,,drop=F]
	DEXSeq_ct(sampleTmp, gffFile, prefix_name, comp, fdr)
}

if (${all}) {
	DEXSeq_ct(sample, gffFile, "${prefix}${mid}.all", 'all', ${fdr})
}

if ("${compare_mode}" == "pairwise" || "${compare_mode}" == "timeseries") {
	if ("${compare_pair}" == "FALSE") {
		compare_data <- as.vector(unique(sample\$conditions))
		#compare_combine <- as.matrix(combn(compare_data, 2))
		#for(i in compare_combine) {
		#	de_twosample(dds, i)
		#}
		len_compare_data <- length(compare_data)
		for(i in 1:(len_compare_data-1)) {
			for(j in (i+1):len_compare_data) {
				tmp_compare <- as.data.frame(
					cbind(sampA=compare_data[i],
					sampB=compare_data[j]))
				de_twosample(sample, gffFile, ${fdr}, tmp_compare)
			}
		}
	}else {
		compare_data <- read.table("${compare_pair}", sep="\t",
		check.names=F, quote='', com='')
		colnames(compare_data) <- c("sampA", "sampB")
		unused <- by(compare_data, 1:nrow(compare_data), function (x)
		de_twosample(sample, gffFile, ${fdr}, x))
	}	
}

file = paste0("${prefix}${mid}.RData")
save.image(file=file)

END

if test "${execute}" == "TRUE";
then
	/bin/rm -f ${all_de}
	R --vanilla < ${prefix}${mid}.r
	touch ${all_de}
fi
