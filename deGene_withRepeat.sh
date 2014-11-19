#!/bin/bash
#set -x
set -e
set -u
#########################################

usage()
{
cat <<EOF
${txtcyn}
Usage (the least parameter):

$0 -f data_file ${txtrst}

${bldblu}Function${txtrst}:

This script is used to detect DE genes by t-test.

${txtbld}OPTIONS${txtrst}:
	-f	Gene expression matrix file.
		${bldred}[Necessary]${txtrst}
	-t	The repeat time of treated sample.[integer, default 3]
	-c	The repeat time of control sample.[integer, default 3]
	-p	The accepted maximum p-value.[default 0.05]
	-d	The accepted maximum fdr.[defaultault 0.3]
	-o	The accepted minimum fold change.
		[log2 based, default 1 means 2 times fold change.]
	-l	Is the data log2 transformed.
		[default TRUE, meaning expression data is already log2
		transformed.
		When FALSE is given, expression data will be transformed
		innerly.]
	-i	If error happends when loading needed packages, plaease give
		TRUE to -i to install all needed ones.
		${bldred}Default [FALSE]${txtrst}
	-r	Run the script[default] or only produce the script[FALSE].
EOF
}

file=
install='FALSE'
run='TRUE'
controlR=3
treatR=3
pvalue=0.05
fdr=0.3
foldc=1
log2="TRUE"
checkNames="TRUE"


while getopts "hf:t:c:p:d:o:l:i:r:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		i)
			install=$OPTARG
			;;
		r)
			run=$OPTARG
			;;
		t)
			treatR=$OPTARG
			;;
		c)
			controlR=$OPTARG
			;;
		p)
			pvalue=$OPTARG
			;;
		d)
			fdr=$OPTARG
			;;
		o)
			foldc=$OPTARG
			;;
		l)
			log2=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done
if [ -z $file ] ; then
	usage
	echo "No -f supplied"
	exit 1
fi

midname="de"

cat <<EOF >$file.${midname}.r
if ($install){
	source("http://bioconductor.org/biocLite.R")
	biocLite(c("genefilter"))
}

library('genefilter')

esetF <- read.table(file="${file}", sep="\t", header=TRUE,
row.names=1, check.names=${checkNames})

esetF <- as.matrix(esetF)

if (! ${log2}){
	esetF <- log2(esetF)
}

Ttest <- rowttests(esetF, as.factor(c(rep(1,$controlR),rep(2,$treatR))))
print('method "BH" gives the false discovery rate ?p.adjust. We declare a collection of 100 genes with a maximum FDR of 0.10 to be differentially expressed (DE), then we expect a maximum of 10 genes to be false positives.')
p.adjust <- p.adjust(Ttest\$p.value, method="BH")
TtestAdj <- cbind(Ttest, p.adjust)
esetFF <- cbind(esetF, TtestAdj)
print("Output all gene expression with t-test results")
write.table(esetFF,
	file="${file}${midname}_${foldc}_${pvalue}_${fdr}.expr.ttest",
	sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
diffExpr <- subset(TtestAdj, abs(TtestAdj\$dm)>=$foldc)
diffExpr <- subset(diffExpr, diffExpr\$p.value<=$pvalue)
diffExpr <- subset(diffExpr, diffExpr\$p.adjust<=$fdr)
write.table(diffExpr, file="${file}${midname}_${foldc}_${pvalue}_${fdr}.deexpr",
	sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
EOF

if [ "$run" = 'TRUE' ];then
	Rscript $file.${midname}.r 
	/bin/rm -f $file.${midname}.r 
fi


