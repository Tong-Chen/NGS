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
	-s	The statistical method you prefer.
		[Default t.test, accept wilcox.test.]	
	-p	The accepted maximum p-value.[default 0.05]
	-d	The accepted maximum fdr.[defaultault 0.3]
	-o	The accepted minimum fold change.
		[log2 based, default 1 means 2 times fold change.]
	-O	Only select DE genes, no statistical needed.
		[Default FALSE, accept TRUE]
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
test_m='t.test'
controlR=3
treatR=3
pvalue=0.05
fdr=0.3
foldc=1
log2="TRUE"
checkNames="TRUE"
onlySelect='FALSE'

while getopts "hf:t:c:s:p:d:o:O:l:i:r:" OPTION
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
		s)
			test_m=$OPTARG
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
		O)
			onlySelect=$OPTARG
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

midname=".de.${test_m}"

cat <<EOF >$file${midname}.r

if (! ${onlySelect}){

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

	if ("${test_m}" == 't.test'){
		sta_test <- rowttests(esetF, as.factor(c(rep(1,$controlR),rep(2,$treatR))))
		sta_test <- sta_test[, 2-3]
		colnames(sta_test) <- c('log2FC', 'p.value')
	} else if ("${test_m}" == 'wilcox.test') {
		rowwilcox.test <- function(x, controlR, treatR) {
			p_fc_list <- c()
			inner_wilcox.test <- function (x, p_fc_list, controlR,treatR){
				control <- x[1:controlR]
				treat <- x[(controlR+1):(controlR+treatR)]
				p <- wilcox.test(control, treat)\$p.value
				fc <- log2(mean(control)/mean(treat))
				if(length(p_fc_list)==0){
					p_fc_list <<- c(fc, p)
				}else {
					p_fc_list <<- rbind(p_fc_list, c(fc, p))
				}
			}
			apply(x, 1, function(x) inner_wilcox.test(x,p_fc_list,controlR,treatR))
			rownames(p_fc_list) <- rownames(x)
			return (as.data.frame(p_fc_list))
		}
		sta_test <- rowwilcox.test(esetF, ${controlR}, ${treatR})
		colnames(sta_test) <- c('log2FC', 'p.value')
	}
	print('method "BH" gives the false discovery rate ?p.adjust. We declare a collection of 100 genes with a maximum FDR of 0.10 to be differentially expressed (DE), then we expect a maximum of 10 genes to be false positives.')
	p.adjust <- p.adjust(sta_test\$p.value, method="BH")
	sta_testAdj <- cbind(sta_test, p.adjust)
	esetFF <- cbind(esetF, sta_testAdj)

	esetFF <- esetFF[order(esetFF\$p.value), ]

	print("Output all gene expression with t-test results")
	write.table(esetFF,
		file="${file}${midname}.expr",
		sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
} else {
	esetFF <- read.table(file="${file}${midname}.expr", sep="\t",
		header=T, row.names=1)
} 

data_len <- ${controlR} + ${treatR}
diffExpr <- subset(esetFF, esetFF\$p.value<=$pvalue)
diffExpr <- subset(diffExpr, diffExpr\$p.adjust<=$fdr)

diffExpr1 <- subset(diffExpr, abs(diffExpr\$log2FC)>=$foldc)
write.table(diffExpr1, file="${file}${midname}_${foldc}_${pvalue}_${fdr}.de.expr",
	sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(diffExpr1[,1:data_len], 
	file="${file}${midname}_${foldc}_${pvalue}_${fdr}.de.expronly",
	sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

diffExpr_up <- subset(diffExpr, diffExpr\$log2FC>=$foldc)
write.table(diffExpr_up,
	file="${file}${midname}_${foldc}_${pvalue}_${fdr}.de.up.expr",
	sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(diffExpr_up[,1:data_len], 
	file="${file}${midname}_${foldc}_${pvalue}_${fdr}.de.up.expronly",
	sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

diffExpr_dw <- subset(diffExpr, diffExpr\$log2FC<=(-1) * $foldc)
write.table(diffExpr_dw,
	file="${file}${midname}_${foldc}_${pvalue}_${fdr}.de.dw.expr",
	sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(diffExpr_dw[,1:data_len], 
	file="${file}${midname}_${foldc}_${pvalue}_${fdr}.de.dw.expronly",
	sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

EOF

if [ "$run" = 'TRUE' ];then
	Rscript $file${midname}.r 
	/bin/rm -f $file${midname}.r 
	sed -i '1 s/^/Gene\t/' ${file}${midname}.expr
	sed -i '1 s/^/Gene\t/' ${file}${midname}_${foldc}_${pvalue}_${fdr}.de.expr
	sed -i '1 s/^/Gene\t/' ${file}${midname}_${foldc}_${pvalue}_${fdr}.de.dw.expr
	sed -i '1 s/^/Gene\t/' ${file}${midname}_${foldc}_${pvalue}_${fdr}.de.up.expr
	sed -i '1 s/^/Gene\t/' ${file}${midname}_${foldc}_${pvalue}_${fdr}.de.expronly
	sed -i '1 s/^/Gene\t/' ${file}${midname}_${foldc}_${pvalue}_${fdr}.de.dw.expronly
	sed -i '1 s/^/Gene\t/' ${file}${midname}_${foldc}_${pvalue}_${fdr}.de.up.expronly
fi


